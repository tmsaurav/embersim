"""

THIS CODE CREATES A MESH WITH MXN OF STRUCTURES WITH REFINEMENT REGION AROUND THEM.
USER CAN SPECIFY A FOREST CANOPY AND WAKE.

TANVIR M. SAURAV

BUSHFIRE RESEARCH GROUP
UNSW CANBERRA

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


gmsh.model.add("valcan") # unique name for mesh
mesh_name = "embersim.msh"
write_mesh = True


# cube parameters

m = 2 # no. of rows of buildings
n = 2 # no. of columns of buildings

clx = 10.0 # length of 'cube'
cly = 10.0 # width of 'cube'
clz = 10.0 # height of 'cube'

cnx = 4     + 1     # cube x-dir elements + 1
cny = 4     + 1     # cube y-dir elements + 1
cnz = 10     + 1     # cube height elements + 1
cnout = 3   # elements out from cube in expansion surfaces
stretch_wall = 1.3;   # stretching near ground
stretch_exp   = 1.3;   # expansion from cube surface


# domain parameters:
    
canopy = 200.0 # length of canopy in the x-direction
apz = 100.0 # length of asset protection zone in the x-direction
lz = 60.0 # height of entire domain
wakelx = 20.0 # wake length

nx = [25] # no. of macro elements in the canopy and apz - multiple numbers denote numel in each layer
nz = [4] # no. of macro elements in the z-dir above cubes - multiple numbers denote numel in each layer

wakenx = [2]; # # no. of macro elements in the wake - multiple numbers denote numel in each layer


#############################################################################
#---------------------------# CUSTOM PARAMETERS  #--------------------------#
#############################################################################

refOffsetx = clx/2.0 # size of refinement zone in x-dir
refOffsety = cly/2.0 # size of refinement zone in y-dir
refOffsetz = 5.0 #clz/2.0 # size of refinement zone in z-dir


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





#############################################################################
#-----------------------# MESH ENGINE - DO NOT TOUCH #----------------------#
#############################################################################

# # ALLOCATE EMPTY LISTS # #

cubepts = []
refpts = []

cubelins = []
reflins = []
creflins = []

cloops = []
refloops = []
intloops = []

csurfs = []
refsurfs = []
intsurfs = []

surfloops = []

vols = []

percubesr = []
percubesl = []

cubefrsurfs = []
cubebksurfs = []
cubeprsurfs = []
cubeplsurfs = []
cubetpsurfs = []

walls = []
inlet = []
outlet = []
topsym = []
periodicr = []
periodicl = []

canvols = []
wakevols = []
cantopvols = []
waketopvols = []
cubetopvols = []

cantopsurfs = []
waketopsurfs = []

fluidvols = []


lx = canopy + apz
lc = 1 

# cube bottom
cubepts.append(gmsh.model.geo.addPoint(0, refOffsety, 0, lc))
cubepts.append(gmsh.model.geo.addPoint(clx, refOffsety, 0, lc))
cubepts.append(gmsh.model.geo.addPoint(clx, refOffsety+cly, 0, lc))
cubepts.append(gmsh.model.geo.addPoint(0, refOffsety+cly, 0, lc))

# cube top
cubepts.append(gmsh.model.geo.addPoint(0, refOffsety, clz, lc))
cubepts.append(gmsh.model.geo.addPoint(clx, refOffsety, clz, lc))
cubepts.append(gmsh.model.geo.addPoint(clx, refOffsety+cly, clz, lc))
cubepts.append(gmsh.model.geo.addPoint(0, refOffsety+cly, clz, lc))


# # refinement bottom
refpts.append(gmsh.model.geo.addPoint(-refOffsetx, 0, 0, lc))
refpts.append(gmsh.model.geo.addPoint(refOffsetx+clx, 0, 0, lc))
refpts.append(gmsh.model.geo.addPoint(refOffsetx+clx, 2*refOffsety+cly, 0, lc))
refpts.append(gmsh.model.geo.addPoint(-refOffsetx, 2*refOffsety+cly, 0, lc))

# # refinement top
refpts.append(gmsh.model.geo.addPoint(-refOffsetx, 0, clz+refOffsetz, lc))
refpts.append(gmsh.model.geo.addPoint(refOffsetx+clx, 0, clz+refOffsetz, lc))
refpts.append(gmsh.model.geo.addPoint(refOffsetx+clx, 2*refOffsety+cly, clz+refOffsetz, lc))
refpts.append(gmsh.model.geo.addPoint(-refOffsetx, 2*refOffsety+cly, clz+refOffsetz, lc))


# cube edges:

# bottom    
cubelins.append(gmsh.model.geo.addLine(cubepts[0], cubepts[1]))
cubelins.append(gmsh.model.geo.addLine(cubepts[1], cubepts[2]))
cubelins.append(gmsh.model.geo.addLine(cubepts[2], cubepts[3]))
cubelins.append(gmsh.model.geo.addLine(cubepts[3], cubepts[0]))

# top
cubelins.append(gmsh.model.geo.addLine(cubepts[4], cubepts[5]))
cubelins.append(gmsh.model.geo.addLine(cubepts[5], cubepts[6]))
cubelins.append(gmsh.model.geo.addLine(cubepts[6], cubepts[7]))
cubelins.append(gmsh.model.geo.addLine(cubepts[7], cubepts[4]))

# sides
cubelins.append(gmsh.model.geo.addLine(cubepts[0], cubepts[4]))
cubelins.append(gmsh.model.geo.addLine(cubepts[1], cubepts[5]))
cubelins.append(gmsh.model.geo.addLine(cubepts[2], cubepts[6]))
cubelins.append(gmsh.model.geo.addLine(cubepts[3], cubepts[7]))

xlc = [1,3,5,7]
ylc = [2,4,6,8]
zlc = [9,10,11,12]


# refinement edges:
 
# bottom    
reflins.append(gmsh.model.geo.addLine(refpts[0], refpts[1]))
reflins.append(gmsh.model.geo.addLine(refpts[1], refpts[2]))
reflins.append(gmsh.model.geo.addLine(refpts[2], refpts[3]))
reflins.append(gmsh.model.geo.addLine(refpts[3], refpts[0]))

# top
reflins.append(gmsh.model.geo.addLine(refpts[4], refpts[5]))
reflins.append(gmsh.model.geo.addLine(refpts[5], refpts[6]))
reflins.append(gmsh.model.geo.addLine(refpts[6], refpts[7]))
reflins.append(gmsh.model.geo.addLine(refpts[7], refpts[4]))

# sides
reflins.append(gmsh.model.geo.addLine(refpts[0], refpts[4]))
reflins.append(gmsh.model.geo.addLine(refpts[1], refpts[5]))
reflins.append(gmsh.model.geo.addLine(refpts[2], refpts[6]))
reflins.append(gmsh.model.geo.addLine(refpts[3], refpts[7]))

xrc = [13,15,17,19]
yrc = [14,16,18,20]
zrc = [21,22,23,24]

# cube and refinement connections:
   
creflins.append(gmsh.model.geo.addLine(cubepts[0], refpts[0]))
creflins.append(gmsh.model.geo.addLine(cubepts[1], refpts[1]))
creflins.append(gmsh.model.geo.addLine(cubepts[2], refpts[2]))
creflins.append(gmsh.model.geo.addLine(cubepts[3], refpts[3]))
creflins.append(gmsh.model.geo.addLine(cubepts[4], refpts[4]))
creflins.append(gmsh.model.geo.addLine(cubepts[5], refpts[5]))
creflins.append(gmsh.model.geo.addLine(cubepts[6], refpts[6]))
creflins.append(gmsh.model.geo.addLine(cubepts[7], refpts[7]))

# cube curve loops and surfaces:



cloops.append(0) # no bottom
csurfs.append(0) # "

# # top
cloops.append(gmsh.model.geo.addCurveLoop([cubelins[4], cubelins[5], cubelins[6], cubelins[7]]))
csurfs.append(gmsh.model.geo.addPlaneSurface([cloops[1]])) # wall
  
# # side

cloops.append(gmsh.model.geo.addCurveLoop([cubelins[0], cubelins[9], -cubelins[4], -cubelins[8]])) 
cloops.append(gmsh.model.geo.addCurveLoop([cubelins[1], cubelins[10], -cubelins[5], -cubelins[9]])) 
cloops.append(gmsh.model.geo.addCurveLoop([cubelins[2], cubelins[11], -cubelins[6], -cubelins[10]])) 
cloops.append(gmsh.model.geo.addCurveLoop([cubelins[3], cubelins[8], -cubelins[7], -cubelins[11]])) 

csurfs.append(gmsh.model.geo.addPlaneSurface([cloops[2]]))
csurfs.append(gmsh.model.geo.addPlaneSurface([cloops[3]]))
csurfs.append(gmsh.model.geo.addPlaneSurface([cloops[4]]))
csurfs.append(gmsh.model.geo.addPlaneSurface([cloops[5]]))
 
# refinement curve loops and surfaces

# # bottom has four sections


refloops.append(gmsh.model.geo.addCurveLoop([reflins[0],-creflins[1],-cubelins[0],creflins[0]]))
refloops.append(gmsh.model.geo.addCurveLoop([reflins[1], -creflins[2], -cubelins[1], creflins[1]]))
refloops.append(gmsh.model.geo.addCurveLoop([reflins[2], -creflins[3], -cubelins[2], creflins[2]]))
refloops.append(gmsh.model.geo.addCurveLoop([reflins[3], -creflins[0], -cubelins[3], creflins[3]]))

refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[0]])) # wall
refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[1]])) # "
refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[2]])) # "
refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[3]])) # "

# # top

refloops.append(gmsh.model.geo.addCurveLoop([reflins[4],reflins[5],reflins[6],reflins[7]]))
refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[4]])) # "

# # sides
refloops.append(gmsh.model.geo.addCurveLoop([reflins[0], reflins[9], -reflins[4], -reflins[8]])) 
refloops.append(gmsh.model.geo.addCurveLoop([reflins[1], reflins[10], -reflins[5], -reflins[9]])) 
refloops.append(gmsh.model.geo.addCurveLoop([reflins[2], reflins[11], -reflins[6], -reflins[10]])) 
refloops.append(gmsh.model.geo.addCurveLoop([reflins[3], reflins[8], -reflins[7], -reflins[11]])) 

refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[5]])) # wall
refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[6]])) # "
refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[7]])) # "
refsurfs.append(gmsh.model.geo.addPlaneSurface([refloops[8]])) # "

# # internal


# # vertical

intloops.append(gmsh.model.geo.addCurveLoop([creflins[0], reflins[8], -creflins[4], -cubelins[8]]))
intloops.append(gmsh.model.geo.addCurveLoop([creflins[1], reflins[9], -creflins[5], -cubelins[9]]))
intloops.append(gmsh.model.geo.addCurveLoop([creflins[2], reflins[10], -creflins[6], -cubelins[10]]))
intloops.append(gmsh.model.geo.addCurveLoop([creflins[3], reflins[11], -creflins[7], -cubelins[11]]))

intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[0]]))
intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[1]]))
intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[2]]))
intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[3]]))

# # top ones

intloops.append(gmsh.model.geo.addCurveLoop([cubelins[4], creflins[5], -reflins[4], -creflins[4]]))
intloops.append(gmsh.model.geo.addCurveLoop([cubelins[5], creflins[6], -reflins[5], -creflins[5]]))
intloops.append(gmsh.model.geo.addCurveLoop([cubelins[6], creflins[7], -reflins[6], -creflins[6]]))
intloops.append(gmsh.model.geo.addCurveLoop([cubelins[7], creflins[4], -reflins[7], -creflins[7]]))

intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[4]]))
intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[5]]))
intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[6]]))
intsurfs.append(gmsh.model.geo.addPlaneSurface([intloops[7]]))



# # add volumes


surfloops.append(gmsh.model.geo.addSurfaceLoop([refsurfs[5], intsurfs[1], csurfs[2], intsurfs[0], refsurfs[0], intsurfs[4]]))
surfloops.append(gmsh.model.geo.addSurfaceLoop([refsurfs[6], intsurfs[2], csurfs[3], intsurfs[1], refsurfs[1], intsurfs[5]]))
surfloops.append(gmsh.model.geo.addSurfaceLoop([refsurfs[7], intsurfs[3], csurfs[4], intsurfs[2], refsurfs[2], intsurfs[6]]))
surfloops.append(gmsh.model.geo.addSurfaceLoop([refsurfs[8], intsurfs[0], csurfs[5], intsurfs[3], refsurfs[3], intsurfs[7]]))

# # top
surfloops.append(gmsh.model.geo.addSurfaceLoop([csurfs[1], refsurfs[4], intsurfs[4], intsurfs[5], intsurfs[6], intsurfs[7]]))

vols.append(gmsh.model.geo.addVolume([surfloops[0]]))
vols.append(gmsh.model.geo.addVolume([surfloops[1]]))
vols.append(gmsh.model.geo.addVolume([surfloops[2]]))
vols.append(gmsh.model.geo.addVolume([surfloops[3]]))

vols.append(gmsh.model.geo.addVolume([surfloops[4]]))


# # set transfinite constraints in x-y-z and surfaces

for i in range(len(xlc)):
    
    gmsh.model.geo.mesh.setTransfiniteCurve(xlc[i], cnx)
    gmsh.model.geo.mesh.setTransfiniteCurve(ylc[i], cny)
    gmsh.model.geo.mesh.setTransfiniteCurve(zlc[i], cnz, "Progression", stretch_wall)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(xrc[i], cnx)
    gmsh.model.geo.mesh.setTransfiniteCurve(yrc[i], cny)
    gmsh.model.geo.mesh.setTransfiniteCurve(zrc[i], cnz, "Progression", stretch_wall)


for i in range(len(creflins)):
    gmsh.model.geo.mesh.setTransfiniteCurve(creflins[i], cnout, "Progression", stretch_exp)


allsurfs = [csurfs,refsurfs,intsurfs]

for i in range(len(allsurfs)):
    
    tmp_surf = allsurfs[i]
    
    for j in range(len(tmp_surf)):
        gmsh.model.geo.mesh.setTransfiniteSurface(tmp_surf[j], "Left")
        gmsh.model.geo.mesh.setRecombine(2, tmp_surf[j])

for i in range(len(vols)):
    gmsh.model.geo.mesh.setTransfiniteVolume(vols[i])


allcubes = [prim_vols for prim_vols in vols]

for i in range(0,n):
    for j in range(0,m):
        for k in range(len(vols)):
            if i == 0 and j == 0:
                break
            tmp_dimtag = gmsh.model.geo.copy([(3, vols[k])])
            allcubes.append(tmp_dimtag[0][1])
            gmsh.model.geo.translate(tmp_dimtag, i*20, j*20, 0)

frontcubes = allcubes[:5*m]
backcubes = allcubes[5*(m*(n-1)):] 

for i in range(0,len(allcubes),5*m):
    percubesr[i:5+i] = allcubes[i:5+i] #########################################################################



for i in range(0,len(allcubes),5*m):
    percubesl[i+5*(m-1):i+5*m] = allcubes[i+5*(m-1):i+5*m]
    
gmsh.model.geo.synchronize()

# # get front surfaces of front cubes

for i in range(len(frontcubes)//5):
    
    tmp_vol = frontcubes[3+5*i]
    
    up, down = gmsh.model.getAdjacencies(3,tmp_vol)
    
    cubefrsurfs.append(down[0])


# # get back surfaces of back cubes

for i in range(len(backcubes)//5):
    
    tmp_vol = backcubes[1+5*i]
    
    up, down = gmsh.model.getAdjacencies(3,tmp_vol)
    
    cubebksurfs.append(down[0])
    
# # periodic right surfaces


   
for i in range(len(percubesr)//5):
    
    tmp_vol = percubesr[5*i]
    
    up, down = gmsh.model.getAdjacencies(3,tmp_vol)
    
    cubeprsurfs.append(down[0])

# # periodic left surfaces
 

for i in range(len(percubesl)//5):
    
    tmp_vol = percubesl[2+5*i]
    
    up, down = gmsh.model.getAdjacencies(3,tmp_vol)
    
    cubeplsurfs.append(down[0])
    
# # top surfaces
 

   
for i in range(len(allcubes)//5):
    
    tmp_vol = allcubes[5*i+4]
    
    up, down = gmsh.model.getAdjacencies(3,tmp_vol)
    
    cubetpsurfs.append(down[1])

# # walls


   
for i in range(len(allcubes)):
    
    tmp_vol = allcubes[i]
    
    up, down = gmsh.model.getAdjacencies(3,tmp_vol)
    
    if i%5 == 4:
        
        walls.append(down[0])
    
    else:
        
        walls.append(down[2])
        walls.append(down[4])



for i in range(len(cubefrsurfs)):
    
    canvols.append(gmsh.model.geo.extrude([(2, cubefrsurfs[i])], -lx+refOffsetx, 0, 0, nx, [1], True))
    wakevols.append(gmsh.model.geo.extrude([(2, cubebksurfs[i])], wakelx-refOffsetx, 0, 0, wakenx, [1], True))

                                           
gmsh.model.geo.synchronize()



for i in range(len(canvols)):
    
    inlet.append(canvols[i][0][1])
    outlet.append(wakevols[i][0][1])
 
    walls.append(canvols[i][2][1])
    walls.append(wakevols[i][2][1])
    
    cantopsurfs.append(canvols[i][4][1])
    waketopsurfs.append(wakevols[i][4][1])
    


periodicr.extend(cubeprsurfs)
periodicl.extend(cubeplsurfs)

for i in range(len(cantopsurfs)):
    
    cantopvols.append(gmsh.model.geo.extrude([(2, cantopsurfs[i])], 0, 0, lz-clz-refOffsetz, nz, [1], True))
    waketopvols.append(gmsh.model.geo.extrude([(2, waketopsurfs[i])], 0, 0, lz-clz-refOffsetz, nz, [1], True))
    
for i in range(len(cubetpsurfs)):
    cubetopvols.append(gmsh.model.geo.extrude([(2, cubetpsurfs[i])], 0, 0, lz-clz-refOffsetz, nz, [1], True))

gmsh.model.geo.synchronize()



for i in range(len(cantopvols)):
    
    inlet.append(cantopvols[i][4][1])
    outlet.append(waketopvols[i][4][1])
    
    topsym.append(cantopvols[i][0][1])
    topsym.append(waketopvols[i][0][1])
    
for i in range(len(cubetopvols)):
    
    topsym.append(cubetopvols[i][0][1])


periodicr.append(canvols[0][3][1])    
periodicr.append(cantopvols[0][5][1])
periodicr.append(wakevols[0][5][1])    
periodicr.append(waketopvols[0][3][1])


periodicl.append(canvols[-1][5][1])    
periodicl.append(cantopvols[-1][3][1])
periodicl.append(wakevols[-1][3][1])
periodicl.append(waketopvols[-1][5][1]) 


for i in range(0,len(cubetopvols),m):
    periodicr.append(cubetopvols[i][2][1])
    periodicl.append(cubetopvols[i+m-1][4][1])


gmsh.model.geo.synchronize()   

gmsh.model.addPhysicalGroup(2, inlet, inlet_id, "inlet")
gmsh.model.addPhysicalGroup(2, outlet, outlet_id, "outlet")
gmsh.model.addPhysicalGroup(2, periodicr, per_r_id, "per_r")
gmsh.model.addPhysicalGroup(2, periodicl, per_l_id, "per_l") 
gmsh.model.addPhysicalGroup(2, walls, wall_id, "wall")
gmsh.model.addPhysicalGroup(2, topsym, top_id, "sym")

getallvols = gmsh.model.getEntities(3)

for i in range(len(getallvols)):
    fluidvols.append(getallvols[i][1])
    
gmsh.model.addPhysicalGroup(3, fluidvols, fluid_id, "fluid")


gmsh.model.geo.synchronize() 
gmsh.model.geo.removeAllDuplicates()

# Generate mesh:
gmsh.model.mesh.generate()
gmsh.model.mesh.removeDuplicateNodes()
 
# Write mesh data:
    
if write_mesh == True:
    gmsh.write(mesh_name)

with open('gridspacing.txt','w') as f:
    
    f.write('### Simulation parameters ###' + '\n' + '\n')
    
    f.write('Cube array: ' + str(m) + ' x ' + str(n) + '\n')
    f.write('Cube dim x - y - z: ' + str(clx) + ' x ' + str(cly) + ' x ' + str(clz) + '\n')
    f.write('Refinement offsets x - y - z: ' + str(refOffsetx) + ' x ' + str(refOffsety) + ' x ' + str(refOffsetz) + '\n')
    f.write('Canopy + APZ in x dir: ' + str(canopy+apz) + '\n')
    f.write('Wake: ' + str(wakelx) + '\n')
    f.write('Domain y-dir: ' + str(m*(cly+2*refOffsety)) + '\n')
    f.write('Domain z-dir: ' + str(lz) + '\n' + '\n')
    
    f.write('### Mesh details ###' + '\n' + '\n')
    
    f.write('Cube surface: dx_c_min: ' + str(clx/(cnx-1)) + '\n')
    f.write('Cube surface: dy_c_min: ' +  str(cly/(cny-1)) + '\n')
    f.write('Cube surface: dz_c_min: ' + str(clz*(stretch_wall-1)/(stretch_wall**(cnz-1)-1)) + '\n')
    f.write('Canopy and APZ stream: dx_can_min: ' + str(-(-lx+refOffsetx)/nx[0]) + '\n')
    f.write('Wake stream: dx_wake_min: ' + str(wakelx/wakenx[0]) + '\n')
    f.write('Global span: dy_min: ' + str((cly+2*refOffsety)/(cny-1)) + '\n')
    f.write('Global wall: dz_min: ' + str((clz+refOffsetz)*(stretch_wall-1)/(stretch_wall**(cnz-1)-1)) + '\n')
    f.write('Refinement stream: dx_ref_min: ' + str((refOffsetx)*(stretch_exp-1)/(stretch_exp**(cnout-1)-1)) + '\n')
    f.write('Refinement span: dy_ref_min: ' + str((refOffsety)*(stretch_exp-1)/(stretch_exp**(cnout-1)-1)) + '\n')
    f.write('Refinement top: dztop_ref_min: ' + str((refOffsetz)*(stretch_exp-1)/(stretch_exp**(cnout-1)-1)) + '\n')
    f.write('Global top: dztop_min: ' + str((lz-refOffsetz-clz)/nz[0]) + '\n')


# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()
 
# It finalize the Gmsh API
#gmsh.finalize()
