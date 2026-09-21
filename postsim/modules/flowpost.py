import pymech as pm
from pymech.vtksuite import hexa2vtk
import numpy as np
import pyvista as pv

def read_wrap_write_pvtk(cases, stat_type):
    """
    Reads and wraps multiple Nek5000 datasets into PyVista format and saves them as .vtu files.
    
    Parameters:
    -----------
    cases : list of str
        List of case names to process.
    stat_type : int
        Type of statistic (1 or 2).
    
    Returns:
    --------
    None
    """
    if not isinstance(cases, list):
        raise TypeError("cases must be a list of strings.")
    if not isinstance(stat_type, int) or stat_type not in [1, 2]:
        raise ValueError("stat_type must be an integer (1 or 2).")
    
    stat_tag = 'STAT1_' if stat_type == 1 else 'STAT2_'
    
    for case in cases:
        datafile = pm.readnek(stat_tag + case + '.f00001')
        dataset = hexa2vtk(datafile)
        dataset = pv.wrap(dataset._vtk_obj)
        dataset.save(stat_tag + case + '.vtu')
        print(f"Dataset saved as {stat_tag + case}.vtu")


def figsize_to_window_size(figsize=(6, 4), dpi=300):
    """Convert matplotlib-style figsize (inches) + dpi into PyVista window_size (pixels)."""
    width_in, height_in = figsize
    width_px = int(width_in * dpi)
    height_px = int(height_in * dpi)
    return [width_px, height_px]


def visualize_grid_slices(
    grid, z_trim=None, slice_normal='z', slice_coords=50,
    scalar='vel', vel_component=0, scalar_name_user=None, 
    scalar_norm=1.0, cmap='viridis', clim=None, zero_contour=True,
    normal_view=True, window_size=None, figsize=(5.905, 2.4), dpi=300,
    show_scalar_bar=False
):
    """
    Visualize slices of an unstructured grid using PyVista after clipping.

    Parameters
    ----------
    grid : pyvista.UnstructuredGrid
        The input grid to visualize.
    z_trim : float, optional
        If set, clips the grid at this z-height before slicing.
    slice_normal : str or tuple
        Direction to slice: 'x', 'y', 'z' or a vector [x, y, z].
    slice_coords : list, optional
        List of coordinates where to create slices. If None, creates 3 evenly spaced slices.
    scalar : str, optional
        Choose between 'vel' or 'temp' for coloring the slices.
    vel_component : int, optional
        If 'vel' is chosen, select component (0: x, 1: y, 2: z).
    scalar_name_user : str, optional
        User-specified scalar name for the colorbar.
    scalar_norm : float, optional
        Normalization factor for scalars.
    cmap : str, optional
        Colormap.
    clim : tuple of (float, float), optional
        Color limits. If None, inferred from data.
    zero_contour : bool, optional
        If True, draws a white contour line at zero.
    normal_view : bool, optional
        If True, sets camera view normal to slice.
    window_size : list, optional
        Pixel resolution of window. If None, computed from figsize and dpi.
    figsize : tuple, optional
        Matplotlib-style figsize in inches.
    dpi : int, optional
        Resolution in dots per inch.
    show_scalar_bar : bool, optional
        Whether to draw a scalar bar. Default False for minimalist PNGs.
    """
    # Clip to z_trim
    if z_trim is not None:
        grid = grid.clip(origin=(0, 0, z_trim), normal=(0, 0, 1))
    
    # Window size
    if window_size is None:
        window_size = figsize_to_window_size(figsize, dpi)

    p = pv.Plotter()
    p.window_size = window_size
    coords = slice_coords

    # Active scalar
    if scalar == 'vel':
        if 'vel' in grid.point_data:
            scalar_name = f"vel_{vel_component}"
            grid.point_data[scalar_name] = grid.point_data['vel'][:, vel_component]
        else:
            raise ValueError("'vel' not found in grid point data.")
    elif scalar == 'temp':
        if 'temp' in grid.point_data:
            scalar_name = 'temp'
        else:
            raise ValueError("'temp' not found in grid point data.")
    else:
        raise ValueError("Invalid scalar. Choose 'vel' or 'temp'.")

    for coord in coords:
        # Slice orientation
        if slice_normal == 'x':
            slice_point = [coord, 0, 0]; normal = [1, 0, 0]; view_direction = 'yz'
            x_min, x_max = coord, coord
            y_min, y_max = 0.0, 150.0
            z_min, z_max = 0.0, min(22.0, z_trim) if z_trim else 22.0
        elif slice_normal == 'y':
            slice_point = [0, coord, 0]; normal = [0, 1, 0]; view_direction = 'xz'
            x_min, x_max = 0.0, 300.0
            y_min, y_max = coord, coord
            z_min, z_max = 0.0, min(22.0, z_trim) if z_trim else 22.0
        else:  # 'z'
            slice_point = [0, 0, coord]; normal = [0, 0, 1]; view_direction = 'xy'
            x_min, x_max = 0.0, 300.0
            y_min, y_max = 0.0, 150.0
            z_min, z_max = coord, coord

        # Slice mesh
        slice_mesh = grid.slice(normal=normal, origin=slice_point)
        if scalar_name_user is None:
            scalar_name_user = scalar_name

        scalars = slice_mesh.point_data.get(scalar_name).astype(float)

        # Normalize
        if scalar_norm is not None and scalar_norm != 0:
            scalars = scalars / scalar_norm

        scalar_range = (scalars.min(), scalars.max()) if scalar_norm else slice_mesh.get_data_range(scalar_name)
        clim = clim if clim is not None else scalar_range

        # Add mesh (no scalar bar by default)
        p.add_mesh(
            slice_mesh, cmap=cmap, scalars=scalars, clim=clim,
            show_scalar_bar=show_scalar_bar,
            scalar_bar_args={
                "title": scalar_name_user,
                "position_x": 0.45,
                "position_y": 0.25,
                "width": 0.15,
                "height": 0.05,
                "vertical": False
            } if show_scalar_bar else {}
        )

        # Zero contour
        if zero_contour:
            contour = slice_mesh.contour(isosurfaces=[0.0])
            if contour.n_points > 0:
                p.add_mesh(contour, color='white', line_width=1.5)

        # Bounds
        slice_bounds = slice_mesh.bounds

        bbox = pv.Box(bounds=slice_bounds)
        p.add_mesh(bbox, color="black", style="wireframe", line_width=1.5)
        
        # Wireframe rectangle
        rectangle = pv.Box(bounds=(x_min, x_max, y_min, y_max, z_min, z_max))
        p.add_mesh(rectangle, color='k', style='wireframe', line_width=1.5)

        # Camera
        if normal_view:
            p.camera_position = view_direction

    return p, scalar_range





def visualize_streamlines(grid, z_trim=None, seed_type='plane', seed_origin=[150, 75, 5], seed_size=[50, 50, 0], 
                          num_seeds=100, integration_time=5.0, max_steps=500, scalar='vel', vel_component=0, 
                          scalar_norm=1.0, cmap='viridis', clim=None, window_size=[1600, 1200]):
    """
    Visualize streamlines in an unstructured grid using PyVista after clipping the data to a height `z_trim`.

    Parameters:
    -----------
    grid : pyvista.UnstructuredGrid
        The input grid to visualize.
    z_trim : float, optional
        If set, clips the grid at this z-height before visualizing.
    seed_type : str, optional
        Type of seed placement ('plane' or 'sphere').
    seed_origin : list, optional
        The center of the seeding region [x, y, z].
    seed_size : list, optional
        Size of the seeding plane [x_size, y_size, z_size].
    num_seeds : int, optional
        Number of seeds for generating streamlines.
    integration_time : float, optional
        Maximum integration time for the streamlines.
    max_steps : int, optional
        Maximum number of steps for streamline integration.
    scalar : str, optional
        Choose between 'vel' or 'temp' for coloring the streamlines.
    vel_component : int, optional
        If 'vel' is chosen, select component (0: x, 1: y, 2: z).
    scalar_norm : float, optional
        User-specified scalar normalization factor.
    cmap : str, optional
        Colormap to use for visualization.
    clim : tuple of (float, float), optional
        User-specified colormap limits.
    window_size : [float, float], optional
        Sets figure size.
    
    Returns:
    --------
    p : pyvista.Plotter
        The PyVista plotter object.
    scalar_range : tuple
        The range of scalar values used in the colormap.
    """

    # Clip the grid if z_trim is specified
    if z_trim is not None:
        grid = grid.clip(origin=(0, 0, z_trim), normal=(0, 0, 1))

    # Create plotter
    p = pv.Plotter()
    p.window_size = window_size

    # Choose active scalar
    if scalar == 'vel':
        if 'vel' in grid.point_data:
            scalar_name = f"vel_{vel_component}"
            grid.point_data[scalar_name] = grid.point_data['vel'][:, vel_component]
        else:
            raise ValueError("'vel' not found in grid point data.")
    elif scalar == 'temp':
        if 'temp' in grid.point_data:
            scalar_name = 'temp'
        else:
            raise ValueError("'temp' not found in grid point data.")
    else:
        raise ValueError("Invalid scalar. Choose 'vel' or 'temp'.")

    # Generate seed points for streamlines
    if seed_type == 'plane':
        x_min, x_max = seed_origin[0] - seed_size[0]/2, seed_origin[0] + seed_size[0]/2
        y_min, y_max = seed_origin[1] - seed_size[1]/2, seed_origin[1] + seed_size[1]/2
        z_min, z_max = seed_origin[2] - seed_size[2]/2, seed_origin[2] + seed_size[2]/2
        
        seed_x = np.random.uniform(x_min, x_max, num_seeds)
        seed_y = np.random.uniform(y_min, y_max, num_seeds)
        seed_z = np.random.uniform(z_min, z_max, num_seeds)
        seeds = np.column_stack((seed_x, seed_y, seed_z))

    elif seed_type == 'sphere':
        radius = seed_size[0] / 2
        phi = np.random.uniform(0, np.pi, num_seeds)
        theta = np.random.uniform(0, 2*np.pi, num_seeds)
        r = np.random.uniform(0, radius, num_seeds)

        seed_x = seed_origin[0] + r * np.sin(phi) * np.cos(theta)
        seed_y = seed_origin[1] + r * np.sin(phi) * np.sin(theta)
        seed_z = seed_origin[2] + r * np.cos(phi)
        seeds = np.column_stack((seed_x, seed_y, seed_z))

    else:
        raise ValueError("Invalid seed type. Choose 'plane' or 'sphere'.")

    # Convert seeds to PyVista polydata
    seed_mesh = pv.PolyData(seeds)

    # Generate streamlines
    streamlines = grid.streamlines_from_source(seed_mesh, integration_direction="both", 
                                               max_time=integration_time, max_steps=max_steps)

    # Get streamline scalars
    scalars = streamlines.point_data.get(scalar_name)

    # Apply user-specified normalization
    if scalar_norm is not None and scalar_norm != 0:
        scalars = scalars / scalar_norm

    # Get color range
    scalar_range = (scalars.min(), scalars.max()) if scalar_norm else streamlines.get_data_range(scalar_name)
    clim = clim if clim is not None else scalar_range

    # Add streamlines to plot
    p.add_mesh(streamlines, scalars=scalars, cmap=cmap, clim=clim, line_width=2.0,
               scalar_bar_args={
                   "title": scalar,
                   "position_x": 0.85,
                   "position_y": 0.1,
                   "width": 0.15,
                   "height": 0.05,
                   "vertical": True
               })

    # Show seed locations
    p.add_mesh(seed_mesh, color="red", point_size=5.0)

    # Show domain boundaries
    p.show_bounds(location='outer', all_edges=True)

    return p, scalar_range
