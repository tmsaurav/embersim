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

#############################################################################
#------------------------# USER-DEFINED PARAMETERS #------------------------#
#############################################################################

lx = 300.0
ly = 150.0
lz = 100.0

case = "valcan"
write_mesh = True

lc = 10

nx = round(lx/lc) + 1
ny = round(ly/lc) + 1
nz = [1,9] #10
layers = [8.0/lz,1]
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

gmsh.model.add("canopy") # unique name for mesh
mesh_name = case + ".msh"
griddata = "grid_" + case + ".txt"

#############################################################################
#-----------------------# MESH ENGINE - DO NOT TOUCH #----------------------#
#############################################################################
botpts = []
botlins = []
botloop = []
botsurf = []

# domain bottom
botpts.append(gmsh.model.geo.addPoint(0, 0, 0, lc))
botpts.append(gmsh.model.geo.addPoint(lx, 0, 0, lc))
botpts.append(gmsh.model.geo.addPoint(lx, ly, 0, lc))
botpts.append(gmsh.model.geo.addPoint(0, ly, 0, lc))

# bottom lines    
botlins.append(gmsh.model.geo.addLine(botpts[0], botpts[1]))
botlins.append(gmsh.model.geo.addLine(botpts[1], botpts[2]))
botlins.append(gmsh.model.geo.addLine(botpts[2], botpts[3]))
botlins.append(gmsh.model.geo.addLine(botpts[3], botpts[0]))

# bottom line loop and surface

botloop.append(gmsh.model.geo.addCurveLoop([botlins[0], botlins[1], botlins[2], botlins[3]]))
botsurf.append(gmsh.model.geo.addPlaneSurface([botloop[0]])) # wall

gmsh.model.geo.mesh.setTransfiniteCurve(botlins[0], nx)
gmsh.model.geo.mesh.setTransfiniteCurve(botlins[1], ny)
gmsh.model.geo.mesh.setTransfiniteCurve(botlins[2], nx)
gmsh.model.geo.mesh.setTransfiniteCurve(botlins[3], ny)

gmsh.model.geo.mesh.setTransfiniteSurface(botsurf[0], "Left")
gmsh.model.geo.mesh.setRecombine(2, botsurf[0])


mesh3d = gmsh.model.geo.extrude([(2, botsurf[0])], 0, 0, lz, nz, layers, True)

gmsh.model.addPhysicalGroup(2, [mesh3d[5][1]], inlet_id, "inlet")
gmsh.model.addPhysicalGroup(2, [mesh3d[3][1]], outlet_id, "outlet")
gmsh.model.addPhysicalGroup(2, [mesh3d[2][1]], per_r_id, "per_r")
gmsh.model.addPhysicalGroup(2, [mesh3d[4][1]], per_l_id, "per_l") 
gmsh.model.addPhysicalGroup(2, [botsurf[0]], wall_id, "wall")
gmsh.model.addPhysicalGroup(2, [mesh3d[0][1]], top_id, "sym")

gmsh.model.addPhysicalGroup(3, [mesh3d[1][1]], fluid_id, "fluid")

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate()
gmsh.model.mesh.removeDuplicateNodes()

if write_mesh == True:
    gmsh.write(mesh_name)

with open(griddata,'w') as f:
    
    f.write('### Simulation parameters ###' + '\n' + '\n')
    
    f.write('Lx: ' + str(lx) + '\n')
    f.write('Ly: ' + str(ly) + '\n')
    f.write('Lz: ' + str(lz) + '\n')
    f.write('Nx: ' + str(nx - 1) + '\n')
    f.write('Ny: ' + str(ny - 1) + '\n')
    f.write('Nz: ' + str(nz) + '\n')

# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()