# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:12:46 2024

@author: z5379427
"""

# Initialization parameters: Don't Touch!

import gmsh
import sys
gmsh.initialize()
gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2) 
gmsh.option.setNumber("Mesh.ElementOrder", 2)

def is_cube_footprint(ii, jj, n, m):
    if ii % 2 == 1 and jj >= 1 and (jj - 1) % 2 == 0:
        col = (ii - 1) // 2
        row = (jj - 1) // 2
        return 0 <= col < n and 0 <= row < m
    return False


#############################################################################
#------------------------# USER-DEFINED PARAMETERS #------------------------#
#############################################################################

case = "ACTESA"
write_mesh = False


# cube parameters
m = 1
n = 1

clx = 20.0 # length of 'cube'
cly = 75.0 # width of 'cube'
clz = 10.0 # height of 'cube'

nz1 = 2

# domain parameters:
    
canopy = 300.0 # length of canopy in the x-direction
apz = 100.0 # length of asset protection zone in the x-direction
lz = 100.0 # height of entire domain
wakelx = 30.0 # wake length

lx = canopy + apz

nz = 9 # no. of macro elements in the z-dir above cubes - multiple numbers denote numel in each layer
lc = 10.0

delta_x = 15.0 #clx/4.0 # size of refinement zone in x-dir
delta_y = 3*5.0 # cly/2.0 # size of refinement zone in y-dir
yedge = 150.0 - m*cly - (m-1)*yspace

#############################################################################
#------------------# ADVANCED PARAMETERS - AVOID TOUCHING #-----------------#
#############################################################################

# # BC IDs # #

inlet_id = 11
outlet_id = 22

per_r_id = 45
per_l_id = 54

wall_id = 90
top_id = 99

fluid_id = 999

gmsh.model.add(case) # unique name for mesh
mesh_name = case + ".msh"
griddata = "grid_" + case + ".txt"

#############################################################################
#-----------------------# MESH ENGINE - DO NOT TOUCH #----------------------#
#############################################################################


# -----------------------
# --- Create Points -----
# -----------------------

pts = []

xpts = [0.0, lx]
for i in range(n):
    xpts.append(xpts[-1] + clx)
    if i < n - 1:
        xpts.append(xpts[-1] + delta_x)
xpts.append(xpts[-1] + wakelx)

ypts = [0.0, yedge / 2]
for i in range(m):
    ypts.append(ypts[-1] + cly)
    if i < m - 1:
        ypts.append(ypts[-1] + delta_y)
ypts.append(ypts[-1] + yedge / 2)

nxs = len(xpts)
nys = len(ypts)

for j in range(nys):
    for i in range(nxs):
        pts.append(gmsh.model.geo.addPoint(xpts[i], ypts[j], 0, lc))


# ------------------------
# --- Create Lines -------
# ------------------------

hlins = []  # horizontal lines
vlins = []  # vertical lines

# Horizontal
for j in range(nys):
    for i in range(nxs - 1):
        start = j * nxs + i
        hlins.append(gmsh.model.geo.addLine(pts[start], pts[start + 1]))

# Vertical
for j in range(nys - 1):
    for i in range(nxs):
        start = j * nxs + i
        vlins.append(gmsh.model.geo.addLine(pts[start], pts[start + nxs]))


lins = hlins + vlins


# ---------------------------------
# --- Create Surfaces and Tags ----
# ---------------------------------

line_loops = []
surfs = []
cube_footprint_surfs = []

num_hlines = (nxs - 1) * nys  # total horizontal lines
num_vlines = nxs * (nys - 1) # sanity check

for j in range(nys - 1):
    for i in range(nxs - 1):
        # Bottom horizontal line (left to right)
        hline_bot = lins[i + j * (nxs - 1)]

        # Top horizontal line (right to left)
        hline_top = -lins[i + (j + 1) * (nxs - 1)]

        # Right vertical line (bottom to top)
        vline_r = lins[num_hlines + (i + 1) + j * nxs]

        # Left vertical line (top to bottom)
        vline_l = -lins[num_hlines + i + j * nxs]

        # Add curve loop in counter-clockwise order
        loop = gmsh.model.geo.addCurveLoop([hline_bot, vline_r, hline_top, vline_l])
        line_loops.append(loop)
        surf = gmsh.model.geo.addPlaneSurface([loop])
        surfs.append(surf)

        if is_cube_footprint(i, j, n, m):
            cube_footprint_surfs.append(surf)


# ------------------------------------
# --- Remove Cube Footprint Surfaces -
# ------------------------------------

for s in cube_footprint_surfs:
    gmsh.model.geo.remove([(2, s)])
    if s in surfs:
        surfs.remove(s)


# ------------------------------------
# --- First layer extrusion ----------
# ------------------------------------

l1canvols =[]
l1wakevols = []
l1canadjvols = []
l1midvols = []

for i in range(2 * m + 1):
    row_offset = i * (2 * n + 1)

    canvol  = 1 + row_offset
    wakevol = (2 * n + 1) + row_offset

    l1canvols.append(gmsh.model.geo.extrude([(2, canvol)], 0, 0, clz, [nz1], [], True))
    l1wakevols.append(gmsh.model.geo.extrude([(2, wakevol)], 0, 0, clz, [nz1], [], True))

    if i % 2 == 0:  # Only even rows
        for j in range(1, 2 * n, 2):  # 2nd, 4th, 6th, ... (j = 1, 3, 5, ...)
            canadjvol = 1 + j + row_offset
            l1canadjvols.append(gmsh.model.geo.extrude([(2, canadjvol)], 0, 0, clz, [nz1], [], True))

    if n > 1:
        for j in range(2, 2 * n, 2):  # j = 2, 4, 6, ..., 2n - 2
            midsurf = 1 + j + row_offset
            l1midvols.append(gmsh.model.geo.extrude([(2, midsurf)], 0, 0, clz, [nz1], [], True))


# ------------------------------------
# ----- First layer meshing ----------
# ------------------------------------

for s in surfs:
    gmsh.model.geo.mesh.setTransfiniteSurface(s, "Left")
    gmsh.model.geo.mesh.setRecombine(2, s)
    

gmsh.model.geo.synchronize() 
gmsh.model.geo.removeAllDuplicates()


# ------------------------------------
# ----- Cube top definitions ---------
# ------------------------------------

csurfl = []
csurfr = []
csurff = []
csurfb =[]
clinl = []
clinr = []
clinf = []
clinb = []


for i in range(m):
    row = 2 * i + 1  # Cube row index in full domain layout

    for j in range(n):
        idx = i * n + j  # row-major index for canadjvols
        mid_idx_l = row * (n - 1) + (j - 1)  # for left
        mid_idx_r = row * (n - 1) + j        # for right

        # ---------- LEFT face ----------
        if j == 0:
            face = l1canvols[2 * i + 1][3][1]
        else:
            face = l1midvols[mid_idx_l][3][1]
        csurfl.append(face)
        up, down = gmsh.model.getAdjacencies(2, face)
        clinl.append(down[2])

        # ---------- RIGHT face ----------
        if j == n - 1:
            face = l1wakevols[2 * i + 1][5][1]
        else:
            face = l1midvols[mid_idx_r][5][1]
        csurfr.append(face)
        up, down = gmsh.model.getAdjacencies(2, face)
        clinr.append(down[2])

        # ---------- FRONT face ----------
        face = l1canadjvols[idx][4][1]
        csurff.append(face)
        up, down = gmsh.model.getAdjacencies(2, face)
        clinf.append(down[2])

        # ---------- BACK face ----------
        face = l1canadjvols[idx + n][2][1]
        csurfb.append(face)
        up, down = gmsh.model.getAdjacencies(2, face)
        clinb.append(down[2])


# ------------------------------------
# -------- Cube top surfaces ---------
# ------------------------------------

ctoploop = []
csurft = []

for i in range(m):
    for j in range(n):
        idx = i * n + j

        loop = gmsh.model.geo.addCurveLoop([
            clinl[idx],   # left
            clinr[idx],   # right
            clinf[idx],   # front
            clinb[idx]    # back
        ])
        ctoploop.append(loop)

        surface = gmsh.model.geo.addPlaneSurface([loop])
        csurft.append(surface)



for s in csurft:
    gmsh.model.geo.mesh.setTransfiniteSurface(s, "Left")
    gmsh.model.geo.mesh.setRecombine(2, s)
    

# ------------------------------------
# ---- Second layer top surfaces -----
# ------------------------------------

l2cansurft = []
l2wakesurft = []
l2canadjsurft = []
l2midsurft = []


for i in range(len(l1canvols)):
    
    l2cansurft.append(l1canvols[i][0][1])   
    l2wakesurft.append(l1wakevols[i][0][1])

if n > 1:
    for i in range(len(l1midvols)):
        l2midsurft.append(l1midvols[i][0][1])

for i in range(len(l1canadjvols)):    
    l2canadjsurft.append(l1canadjvols[i][0][1])  


l2surfs = l2cansurft + l2canadjsurft + l2midsurft + l2wakesurft + csurft


# ------------------------------------
# ------ Second layer extrusion ------
# ------------------------------------

l2canvols = []
l2wakevols = []
l2canadjvols = []
l2midvols = []

l2cvols = []
l2vols = []

# Extrude cube tops
for s in l2cansurft:
    l2canvols.append(gmsh.model.geo.extrude([(2, s)], 0, 0, lz - clz, [nz], [], True))

# Extrude canonical adjacents (front/back)
for s in l2canadjsurft:
    l2canadjvols.append(gmsh.model.geo.extrude([(2, s)], 0, 0, lz - clz, [nz], [], True))

# Extrude mid cubes (between canonical and wake)
if n > 1:
    for s in l2midsurft:
        l2midvols.append(gmsh.model.geo.extrude([(2, s)], 0, 0, lz - clz, [nz], [], True))

# Extrude wake tops
for s in l2wakesurft:
    l2wakevols.append(gmsh.model.geo.extrude([(2, s)], 0, 0, lz - clz, [nz], [], True))

# Extrude custom top surfaces (csurft from face recombination)
for s in csurft:
    l2cvols.append(gmsh.model.geo.extrude([(2, s)], 0, 0, lz - clz, [nz], [], True))

l2vols = l2canvols + l2canadjvols + l2midvols + l2wakevols + l2cvols

gmsh.model.geo.synchronize() 
gmsh.model.geo.removeAllDuplicates()




# ------------------------------------
# ---------- BC tagging --------------
# ------------------------------------


# ------------------------------------
# -------------- Inlet ---------------
# ------------------------------------

inlet = []

for i in range(len(l1canvols)):
    
    inlet.append(l1canvols[i][5][1])
    inlet.append(l2canvols[i][5][1])


# ------------------------------------
# -------------- Outlet --------------
# ------------------------------------

outlet = []

for i in range(len(l1wakevols)):

    outlet.append(l1wakevols[i][3][1])
    outlet.append(l2wakevols[i][3][1])


# ------------------------------------
# -------------- Top -----------------
# ------------------------------------

top = []

for i in range(len(l2vols)):

    top.append(l2vols[i][0][1])

# ------------------------------------
# ------------ Walls -----------------
# ------------------------------------

wall = []

wall = surfs + csurfl + csurfr + csurff + csurfb + csurft 

# ------------------------------------
# ----------- Periodic ---------------
# ------------------------------------

per_r = [
    l1canvols[0][2][1],
    *( [vol[2][1] for vol in l1midvols[:n - 1]] if n > 1 else [] ),
    l1wakevols[0][2][1],
    *[vol[2][1] for vol in l1canadjvols[:n]],
    l2canvols[0][2][1],
    *( [vol[2][1] for vol in l2midvols[:n - 1]] if n > 1 else [] ),
    l2wakevols[0][2][1],
    *[vol[2][1] for vol in l2canadjvols[:n]]
]

per_l = [
    l1canvols[-1][4][1],
    *( [vol[4][1] for vol in l1midvols[-(n - 1):]] if n > 1 else [] ),
    l1wakevols[-1][4][1],
    *[vol[4][1] for vol in l1canadjvols[-n:]],
    l2canvols[-1][4][1],
    *( [vol[4][1] for vol in l2midvols[-(n - 1):]] if n > 1 else [] ),
    l2wakevols[-1][4][1],
    *[vol[4][1] for vol in l2canadjvols[-n:]]
]


fluidvols = []

getallvols = gmsh.model.getEntities(3)

for i in range(len(getallvols)):
    fluidvols.append(getallvols[i][1])
    
gmsh.model.addPhysicalGroup(2, inlet, inlet_id, "inlet")
gmsh.model.addPhysicalGroup(2, outlet, outlet_id, "outlet")
gmsh.model.addPhysicalGroup(2, per_r, per_r_id, "periodic_r")
gmsh.model.addPhysicalGroup(2, per_l, per_l_id, "periodic_l")
gmsh.model.addPhysicalGroup(2, top, top_id, "top")
gmsh.model.addPhysicalGroup(2, wall, wall_id, "wall")
gmsh.model.addPhysicalGroup(3, fluidvols, fluid_id, "fluid")

gmsh.model.geo.synchronize() 
gmsh.model.geo.removeAllDuplicates()
gmsh.model.mesh.generate()
gmsh.model.mesh.removeDuplicateNodes()


if write_mesh == True:
    gmsh.write(mesh_name)
    
    
# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()

