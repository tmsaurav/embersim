import vtk
from vtk.util.numpy_support import vtk_to_numpy #thats what you need 
import numpy as np
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon, Rectangle
import matplotlib.cm as cm
import matplotlib.colors as colors

def read_vtu_data(file_name):
    
    # Read the source file.
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()

    fldd = output.GetFieldData()
    time = fldd.GetAbstractArray(0).GetComponent(0,0)
    
    potential   = {} 
    flag = 0
    i    = 1
    nfield = 0
    nflag = 0
    while flag == 0:
        nfield += 1
        try:
            potential[str(i)] = vtk_to_numpy(output.GetPointData().GetArray("y0"+str(i)+"   "))    
        except:
            flag = 1
            nflag += 1
        i = i + 1

    flag = 0
    ii  =  1
    while flag == 0:
        nfield += 1
        try:
            potential[str(i)] = vtk_to_numpy(output.GetPointData().GetArray("rpro0"+str(ii)))    
        except:
            flag = 1
            nflag += 1
        i = i + 1
        ii = ii + 1

    flag = 0
    ii  =  1
    while flag == 0:
        nfield += 1
        try:
            potential[str(i)] = vtk_to_numpy(output.GetPointData().GetArray("tag0"+str(ii)+" "))   
        except:
            flag = 1
            nflag += 1
        i = i + 1
        ii = ii + 1
        
    flag = 0 
    for i in potential.keys():
        if flag == 0:
            x       = potential[i]
            flag    = 1
        else:    
            x       = np.vstack((x,potential[i]))

    nfield = nfield - nflag
    return x.transpose(), time, nfield
        

def readsort_particle_data(pathtodir,start,end,npart):
    _, _, nfield = read_vtu_data(pathtodir + "par{:05d}.vtu".format(start))
    
    all_sorted = np.zeros((end-start+1,npart,nfield+1)) # field+1 for time
    jj = 0
    for ii in range(start,end+1):
        file = pathtodir + "par{:05d}.vtu".format(ii)
        data, time, _ = read_vtu_data(file)
        data = data[data[:,-1].argsort()]
        data = data[data[:,-2].argsort(kind='stable')]
        all_sorted[jj,:data.shape[0],:nfield] = data
        all_sorted[jj,:,-1] = time
        jj += 1

    # Extract the particle identifiers (last two columns)
    particle_ids = all_sorted[:, :, -3:-1]

    # Get unique particle identifiers across all timesteps
    unique_particle_ids = np.unique(particle_ids.reshape(-1, 2), axis=0)
    unique_particle_ids = unique_particle_ids[~np.all(unique_particle_ids == [0, 0], axis=1)]

    # Create a mapping from particle identifier tuple to index
    particle_id_to_index = {tuple(pid): idx for idx, pid in enumerate(unique_particle_ids)}

    # Initialize the new array with NaNs
    N_unique_particles = len(unique_particle_ids)
    sorted_data = np.full((all_sorted.shape[0], N_unique_particles, all_sorted.shape[2]), np.nan)

    # Set identifiers in new_array
    for idx, pid in enumerate(unique_particle_ids):
        sorted_data[:, idx, -3:-1] = pid

    # Populate the new array with properties from data array
    for t in range(all_sorted.shape[0]):
        for p in range(all_sorted.shape[1]):
            pid = tuple(all_sorted[t, p, -3:-1])  # Extract the particle identifier tuple
            if pid in particle_id_to_index:
                idx = particle_id_to_index[pid]
            # Set the properties in the new array
                sorted_data[t, idx, :-3] = all_sorted[t, p, :-3]
                sorted_data[t, idx, -1] = all_sorted[t, p, -1]  # Copy the last property as well
        
    return sorted_data


def create_polygon_array(start_x, start_y, rows, cols, clx, cly, xspace, yspace):
    '''
        start_x = 0  # x-coordinate of the first polygon's bottom-left corner
        start_y = 0  # y-coordinate of the first polygon's bottom-left corner
        rows = 3     # number of rows in the array
        cols = 3     # number of columns in the array
        clx = 2        # width of each polygon
        cly = 3        # height of each polygon
        xspace = 1   # horizontal spacing between polygons
        yspace = 1   # vertical spacing between polygons
    '''
    
    polygons = []
    for i in range(rows):
        for j in range(cols):
            bottom_left_x = start_x + j * (clx + xspace)
            bottom_left_y = start_y + i * (cly + yspace)
            
            polygon = Polygon([
                (bottom_left_x, bottom_left_y),
                (bottom_left_x + clx, bottom_left_y),
                (bottom_left_x + clx, bottom_left_y + cly),
                (bottom_left_x, bottom_left_y + cly)
            ])
            
            polygons.append(polygon)
    return polygons


def count_polygon_hits(particle_positions, polygons, z_thresh):
    """
    Identifies particles that are present at timestep t-1 but become invalid (NaN) at timestep t,
    and determines which polygon they were inside before disappearing.
    Tracks cumulative counts of such particle losses per polygon over time.

    Parameters:
    -----------
    particle_positions : np.ndarray
        3D numpy array of shape (time, particles, location), where location contains [x, y, z].
        NaN in x or y indicates the particle is no longer present.

    polygons : list of Shapely Polygon
        List of polygons defining the areas of interest (e.g., cubes).

    z_thresh : float
        Only particles with z <= z_thresh are considered for counting inside polygons.

    Returns:
    --------
    losing_positions : list of tuples
        Each tuple contains (timestep, particle_index, polygon_index) indicating when and where a particle was lost.

    cumulative_counts : np.ndarray
        2D array of shape (time, num_polygons) tracking the cumulative number of lost particles per polygon over time.
    """

    losing_positions = []
    num_polygons = len(polygons)
    cumulative_counts = np.zeros((particle_positions.shape[0], num_polygons), dtype=int)

    # Iterate over timesteps (starting from 1 since we need t-1)
    for t in range(1, particle_positions.shape[0]):
        # Mask: particles that were valid in t-1
        valid_previous = ~np.isnan(particle_positions[t - 1, :, :2]).any(axis=1)
        # Mask: particles that are invalid in t
        valid_current = np.isnan(particle_positions[t, :, :2]).any(axis=1)
        # Mask: particles that disappeared between t-1 and t
        losing_mask = valid_previous & valid_current
        losing_indices = np.where(losing_mask)[0]

        # Check each lost particle
        for index in losing_indices:
            x_prev, y_prev, z_prev = particle_positions[t - 1, index]
            if np.isnan(z_prev):
                continue  # Skip if z is also NaN for safety

            # Check in which polygon the particle was located (with z check)
            for i, polygon in enumerate(polygons):
                if shapely.contains_xy(polygon, x_prev, y_prev) and z_prev <= z_thresh:
                    losing_positions.append((t, index, i, x_prev, y_prev, z_prev))
                    cumulative_counts[t, i] += 1
                    break  # Stop after first matching polygon

        # Carry forward previous cumulative counts to maintain accumulation over time
        if t > 1:
            cumulative_counts[t] += cumulative_counts[t - 1]

    return losing_positions, cumulative_counts


def heightclip_particle_data(sorted_data,low_th,up_th,filtertype):

    if filtertype == 'leq':
        mask = sorted_data[:, :, 2] <= low_th
    elif filtertype == 'geq':
        mask = sorted_data[:, :, 2] >= up_th
    elif filtertype == 'between':
        mask = (sorted_data[:, :, 2] >= low_th) & (sorted_data[:, :, 2] <= up_th)
    else:
        raise ValueError("Choose either 'leq' or 'geq' or 'between' filtertype")
        return

    sorted_data_clipped = np.full_like(sorted_data, np.nan)

    for i in range(sorted_data_clipped.shape[0]):
        sorted_data_clipped[i, mask[i], :] = sorted_data[i, mask[i], :]
    
    return sorted_data_clipped


def remove_periodic_x(data,lower_bound,upper_bound):
    
    # Compute the change in positions between timesteps
    position_diff = np.diff(data, axis=0)
    
    # Find particles with jumps exceeding half the domain (indicative of periodicity)
    domain_width = upper_bound - lower_bound
    wrap_condition = np.abs(position_diff) > (domain_width / 2)
    
    # Identify particles that wrapped at any timestep
    wrapped_particles = np.where(wrap_condition.any(axis=0))[0]
    
    # Create a copy of the data to store non-periodic positions
    data_nonperiodic = np.copy(data)
    
    # Loop through wrapped particles and set subsequent locations to NaN
    for particle in wrapped_particles:
        # Find the first timestep where wrapping occurred
        wrap_timestep = np.where(wrap_condition[:, particle])[0][0] + 1  # +1 due to diff
        # Set all subsequent positions to NaN
        data_nonperiodic[wrap_timestep:, particle] = np.nan

    return data_nonperiodic



def bin_average_particle_data(data, nx, ny,
                                   x_bounds,
                                   y_bounds,
                                   output_file):
    """
    Bin and average particle data in 2D spatial bins over time.

    Parameters:
    -----------
    data : ndarray of shape (T, N, P)
        Time-series particle data with at least 13 properties.
    nx, ny : int
        Number of bins in x and y directions.
    x_bounds, y_bounds : tuple
        Min and max bounds for x and y dimensions.
    output_file : str
        Path to save the compressed .npz output.

    Returns:
    --------
    None. Saves output to `output_file`.
    """

    T, N, P = data.shape
    x_edges = np.linspace(x_bounds[0], x_bounds[1], nx + 1)
    y_edges = np.linspace(y_bounds[0], y_bounds[1], ny + 1)

    counts = np.zeros((T, nx, ny), dtype=int)
    
    field_names = [
        "vx", "vy", "vz", "rho", "dp", "psip",
        "ux", "uy", "uz", "dudz"
    ]
    field_indices = dict(zip(field_names, range(3, 13)))  # columns 3 to 12

    avg_fields = {
        name: np.full((T, nx, ny), np.nan, dtype=np.float32)
        for name in field_names
    }

    for t in range(T):
        x = data[t, :, 0]
        y = data[t, :, 1]

        x_idx = np.digitize(x, x_edges) - 1
        y_idx = np.digitize(y, y_edges) - 1
        valid = (x_idx >= 0) & (x_idx < nx) & (y_idx >= 0) & (y_idx < ny)

        field_data = {
            name: data[t, :, idx] for name, idx in field_indices.items()
        }

        for i in range(N):
            if valid[i]:
                xi, yi = x_idx[i], y_idx[i]
                counts[t, xi, yi] += 1
                for name in field_names:
                    if np.isnan(avg_fields[name][t, xi, yi]):
                        avg_fields[name][t, xi, yi] = field_data[name][i]
                    else:
                        avg_fields[name][t, xi, yi] += field_data[name][i]

        mask = counts[t] > 0
        for name in field_names:
            avg_fields[name][t][mask] /= counts[t][mask]

    # Save all fields
    np.savez_compressed(output_file,
                        x_edges=x_edges,
                        y_edges=y_edges,
                        counts=counts,
                        **{f"avg_{name}": array for name, array in avg_fields.items()})


def init_case_params(case_id):
    config = {
        5:    {'cases': ['I80Y5_r3', 'I80Y5_r4', 'I100Y5_r3', 'I100Y5_r4'], 'rows': 6, 'cols': 2, 'xspace': 15.0, 'yspace': 5.0},
        10:   {'cases': ['I80Y10_r3', 'I80Y10_r4', 'I100Y10_r3', 'I100Y10_r4'], 'rows': 5, 'cols': 2, 'xspace': 15.0, 'yspace': 10.0},
        15:   {'cases': ['I80Y15_r3', 'I80Y15_r4', 'I100Y15_r3', 'I100Y15_r4'], 'rows': 5, 'cols': 2, 'xspace': 15.0, 'yspace': 15.0},
        'norm': {'cases': ['norm_I80Y10_r3', 'norm_I80Y10_r4', 'norm_I100Y10_r3', 'norm_I100Y10_r4'], 'rows': 5, 'cols': 2, 'xspace': 15.0, 'yspace': 10.0},
    }

    if case_id not in config:
        raise ValueError(f"Invalid case_id {case_id}. Valid options: {list(config.keys())}")
    
    conf = config[case_id]
    return conf['cases'], conf['rows'], conf['cols'], conf['xspace'], conf['yspace']



####################################### TEST FUNCTIONS ########################################

def checkThrough(polygons,start,end,npart,all_sc,bias = 0):
    
    # bias shifts the start of time range
    
    hits = np.zeros((end-start+1,len(polygons)))
    
    for npi in range(0,npart):    
        breaker = False
        
        for ii in range(bias,end-start+1):
            coords = all_sc[ii,npi,:2]
            coords = Point(coords)
            #print(str(npi) + '/' + str(npart) + ', ' + str(ii))
            for jj in range(0,len(polygons)):
                cond = (polygons[jj].contains(coords))
                
                if cond is True:
                    print(str(npi) + '/' + str(npart) + ', ' + str(ii) + ', ' + str(jj))
                    all_sc[ii:,npi,:2] = -99999
                    all_sc[ii:,npi,19] = jj  ## hit marker
                    hits[ii,jj] = hits[ii,jj] + 1
                    breaker = True
                    break
            if breaker:
                break
    return all_sc, hits 

def check_in_box(x, y, z, region_bounds):
    """
    Check if points (x, y, z) are within a 3D axis-aligned bounding box.

    Parameters:
        x, y, z : float or np.ndarray
            Coordinates of the particles.
        region_bounds : dict
            Dictionary with keys 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'.

    Returns:
        mask : bool or np.ndarray of bool
            True where the point(s) are inside the region.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    
    return ((region_bounds['xmin'] <= x) & (x <= region_bounds['xmax']) &
            (region_bounds['ymin'] <= y) & (y <= region_bounds['ymax']) &
            (region_bounds['zmin'] <= z) & (z <= region_bounds['zmax']))

def plot_binned_quantity(ax,npz_file, quantity_key, time_mode='mean', tstep=None,
                         cmap='RdBu', vmin=None, vmax=None, label=None, title=None,
                         ):
    data = np.load(npz_file)
    x_edges = data['x_edges']
    y_edges = data['y_edges']
    quantity = data[quantity_key]

    x_min, x_max = x_edges[0], x_edges[-1]
    y_min, y_max = y_edges[0], y_edges[-1]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = ax.figure
        
    if time_mode == 'slice':
        assert tstep is not None, "tstep must be specified for 'slice' mode"
        plot_data = quantity[tstep]
    elif time_mode == 'mean':
        plot_data = np.nanmean(quantity, axis=0)
    else:
        raise ValueError("Invalid time_mode. Use 'mean' or 'slice'.")

    im = ax.imshow(
        plot_data.T, origin='lower',
        extent=[x_min, x_max, y_min, y_max],
        aspect='equal', cmap=cmap,
        vmin=vmin, vmax=vmax
    )
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(title if title else f'{time_mode.capitalize()} of {quantity_key}')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(label if label else quantity_key)

    plt.tight_layout()
    
    return fig, ax, im, plot_data



def plot_domain_with_cubes(ax, polygs, rows, cols,
                           cube_facecolors=None,
                           canopy_facecolor=(0, 0.8, 0, 0.2),
                           show_grid=True,
                           zorder=10):  # 👈 Add zorder as a parameter
    
    canopy = Rectangle((0, 0), 300, 150,
                       linewidth=0.8,
                       edgecolor='w',
                       facecolor='none',
                       linestyle = '--',
                       zorder=zorder-4)  # 👈 Same for canopy
    #ax.add_patch(canopy)
    
    for i in range(rows * cols):
        if cube_facecolors is None:
            facecolor = (1, 1, 1, 1)
        elif isinstance(cube_facecolors[0], (float, int)):
            facecolor = cube_facecolors
        else:
            facecolor = cube_facecolors[i]

        cube_patch = MplPolygon(polygs[i].exterior.coords,
                                facecolor=facecolor,
                                edgecolor='k' if show_grid else 'none',
                                linewidth=1,
                                zorder=zorder)  # 👈 Force cubes on top
        ax.add_patch(cube_patch)




def cube_facecolors_from_counts(counts_at_t, cmap='Reds', vmin=None, vmax=None, log_scale=False):
    if log_scale:
        norm = colors.LogNorm(vmin=max(vmin, 1e-2), vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin if vmin is not None else 0,
                                vmax=vmax if vmax is not None else np.max(counts_at_t))
    scalar_cmap = cm.get_cmap(cmap)
    cube_facecolors = [scalar_cmap(norm(count)) for count in counts_at_t]
    return cube_facecolors

