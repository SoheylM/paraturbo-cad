#adapted from May 11 version of takingrotor.py

import pickle
import numpy as np
import cadquery as cq
import os
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions
from cq_warehouse import *
from cadquery import exporters

from cqmore import Workplane
from cqmore.polyhedron import gridSurface

#open file
file = open('Element_23_06_05_v2.pickle', 'rb')

# dump info to that file
Element = pickle.load(file)

#close file
file.close()



file1 = open('Element_23_06_05_v6.pickle', 'rb')

# dump info to that file
Element1 = pickle.load(file1)

#close file
file1.close()




cwf = os.getcwd().replace("\\", "/")

# Use this in CQ-Editor
Rotor = cq.importers.importStep(cwf + '/Rotor.stp').translate((-0.5,1,-15))

# Use this in VS Code
# Rotor = cq.importers.importStep(cwf + '/HGJB/Rotor.stp')


#Rotor = Rotor.rotate((0,0,0),(1,0,0),270)

Laenge = Element['Laenge']
#change list into numpy array
Laenge = np.array(Laenge)
#change units from m to mm to avoid hollow visual effect

Laenge = 1000*Laenge 
#multiplying Laenge streeeetches the part out, so don't
DI1 = 1000*np.array(Element['DI1'])
DI2 = 1000*np.array(Element['DI2'])
DI3 = 1000*np.array(Element['DI3'])
DA1 = 1000*np.array(Element['DA1'])
DA2 = 1000*np.array(Element['DA2'])
DA3 = 1000*np.array(Element['DA3'])
#print('type DI1', type(DI1))

#find position of HGJBs
sys_pos = Element['sys_pos']
pos_hgjb1 = sys_pos['pos_hgjb1']
pos_hgjb2 =sys_pos['pos_hgjb2']


#Begin code adapted from Christophe's Matlab
#Parameters 
N_HG = 28 #number of grooves generally between 26 - 30
alpha_HG = Element['parameters']['hgjb1']['alpha']#0.68 #given
beta_HG = Element['parameters']['hgjb1']['beta'] #-135 #given-135 #given 
beta_HG = beta_HG*pi/180
gamma_HG = Element['parameters']['hgjb1']['gamma'] #0.89 #given
h_gr = Element['parameters']['hgjb1']['hg'] #16 #groove depth given in micrometers
h_rr = 9 #clearance on radiu given in micrometers
D = 16 #on drawing [mm]
L = Laenge[pos_hgjb1]#28 #length of HGJB on drawing [mm]
L_land=L-(gamma_HG*L) #Value for CAD
L=L+0.8 #oversized length for safety generally between 0.6 - 1
Spiral_step = pi*D*tan(beta_HG)
Spiral_height = L/2
a_HG = (pi*D*alpha_HG)/N_HG #mm
a_HG_plus_b_HG = a_HG/alpha_HG #mm
h_rr_tot = h_rr*2 #diametral clearance given in micrometers
#end code adapted from Christophe's Matlab

# number of parallelograms discretized
n_parall = 10

# percentage epsilon of length to extend to avoid surfaces between parallelograms
eps_perc = 0.005

#find distance to first HGJB
dist1 = 0
for d in range(pos_hgjb1):
    dist1=dist1+Laenge[d]
    
#find distance to first HGJB
dist2 = 0
for d in range(pos_hgjb2):
    dist2=dist2+Laenge[d]
    
#find distance to center of first HGJB
DistCenter1=dist1+L/2

#find distance to center of second HGJB
DistCenter2=dist2+L/2

#define separation angle between grooves
sepang=360/N_HG

#length between parallelogram verticals
LenBetwVert = (L/2 - L_land/2)/n_parall

#Angle between parallelogram diagonal and reference horizontal
Betaprime = abs(beta_HG)- pi/2# - 90

#gap between horizontal leaving from parallelogram lowest 
#corner and parallelogram vertical

gap = LenBetwVert*tan(Betaprime)
gap_spiral = (gap*Spiral_step)/LenBetwVert

#create cylinder parameters
CylLen1 = Laenge[pos_hgjb1]
CylRadOut1= DA3[pos_hgjb1]/2 #1.01*

hgjb1_x_first = Element['parameters']['hgjb1']['x_first_curve']
hgjb1_y_first = Element['parameters']['hgjb1']['y_first_curve']
hgjb1_z_first = Element['parameters']['hgjb1']['z_first_curve']
coords_hgjb1_first = list(zip(hgjb1_x_first, hgjb1_y_first, hgjb1_z_first))


# Adding the first point of the curve to the end of its point array to close it
# coords_hgjb1_first.append(coords_hgjb1_first[0])














from cqmore import Workplane
from cqmore.polyhedron import gridSurface

inner_1 = Element1['parameters']['hgjb1']['x_first_curve'][0:10]
inner_2 = Element1['parameters']['hgjb1']['y_first_curve'][0:10]
inner_3 = Element1['parameters']['hgjb1']['z_first_curve'][0:10]
inner1 = list(zip(inner_1,inner_2,inner_3))

inner_11 = Element1['parameters']['hgjb1']['x_first_curve'][10:30]
inner_21 = Element1['parameters']['hgjb1']['y_first_curve'][10:30]
inner_31 = Element1['parameters']['hgjb1']['z_first_curve'][10:30]
inner11 = list(zip(inner_11,inner_21,inner_31))

inner11.reverse()

inner = [inner1,inner11]

surface = Workplane().polyhedron(*gridSurface(inner,1))

show_object(surface)

Rotor = Rotor.cut(surface)

show_object(Rotor)
















# inner_111 = Element1['parameters']['hgjb1']['x_first_surface'][101:150]
# inner_211 = Element1['parameters']['hgjb1']['y_first_surface'][101:150]
# inner_311 = Element1['parameters']['hgjb1']['z_first_surface'][101:150]
# inner111 = list(zip(inner_111,inner_211,inner_311))



# inner = [coords_hgjb1_first]

hgjb1_x_second = Element['parameters']['hgjb1']['x_second_curve']
hgjb1_y_second = Element['parameters']['hgjb1']['y_second_curve']
hgjb1_z_second = Element['parameters']['hgjb1']['z_second_curve']
coords_hgjb1_second = list(zip(hgjb1_x_second, hgjb1_y_second, hgjb1_z_second))

hgjb2_x_first = Element['parameters']['hgjb2']['x_first_curve']
hgjb2_y_first = Element['parameters']['hgjb2']['y_first_curve']
hgjb2_z_first = Element['parameters']['hgjb2']['z_first_curve']
coords_hgjb2_first = list(zip(hgjb2_x_first, hgjb2_y_first, hgjb2_z_first))

hgjb2_x_second = Element['parameters']['hgjb2']['x_second_curve']
hgjb2_y_second = Element['parameters']['hgjb2']['y_second_curve']
hgjb2_z_second = Element['parameters']['hgjb2']['z_second_curve']
coords_hgjb2_second = list(zip(hgjb2_x_second, hgjb2_y_second, hgjb2_z_second))




###################################################



# surface = Workplane().splineApproxSurface([cq.Vector(p) for p in inner],0.2)



# def paraboloid(x, y):
#     return (x, y, ((y ** 2) - (x ** 2)) / 4)

# min_value = -30
# max_value = 30
# step = 5
# thickness = 0.1

# points = [[
#         paraboloid(x / 10, y / 10) 
#     for y in range(min_value, max_value + step, step)
# ] for x in range(min_value, max_value + step, step)]

# surface = Workplane().splineApproxSurface(points, thickness)

# surface = Workplane().splineApproxSurface(inner,0.1)



###################################################







# edge_hgjb1_first_1 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][0:30])
# edge_hgjb1_first_2 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][50:60])
# edge_hgjb1_first_5 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][30:60])

# coords_hgjb1_first.reverse()


# edge_hgjb1_first_3 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][50:80])
# edge_hgjb1_first_4 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][150:200])
# edge_hgjb1_first_6 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][80:100])



# edge_hgjb1_first_2 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][50:60])




# edge_hgjb1_first_4 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][150:200])







###################################################
# edge_hgjb1_second_1 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_second][0:30])
# coords_hgjb1_second.reverse()
# edge_hgjb1_second_3 = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_second][50:80])

###################################################


# .close()


# coordinates = []
# coordinates.extend(coords_hgjb1_first[0:50])
# coordinates.extend(coords_hgjb1_first[51:100])
# coordinates.extend(coords_hgjb1_first[101:150])
# coordinates.extend(coords_hgjb1_first[151:200])


# show_object(edge_hgjb1_first_1)
# show_object(edge_hgjb1_first_2)
# show_object(edge_hgjb1_first_3)
# show_object(edge_hgjb1_first_4)


# edge_hgjb1_first = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_first][0:200]).close().clean()
# edge_hgjb1_second = cq.Edge.makeSplineApprox([cq.Vector(p) for p in coords_hgjb1_second][0:200]).close().clean()

# show_object(edge_hgjb1_first)

# result1 = cq.Workplane("front").polyline(coords_hgjb1_first[0:50])
# show_object(result1)
# result2 = cq.Workplane("front").polyline(coords_hgjb1_first[50:100])
# show_object(result2)
# result3 = cq.Workplane("front").polyline(coords_hgjb1_first[100:150])
# show_object(result3)
# result4 = cq.Workplane("front").polyline(coords_hgjb1_first[150:200])
# show_object(result4)





# result5 = cq.Face.makeNSidedSurface([result1,result2,result3,result4],[])

# face_hgjb1_first = cq.Wire.makeNonPlanarFace([edge_hgjb1_first],surfacePoints=[cq.Vector(p) for p in inner])
# show_object(face_hgjb1_first)

# planar_grid = cq.Wire.makePolygon([cq.Vector(v) for v in inner])
# show_object(planar_grid)

# coordinates1 = []
# coordinates1.extend(coords_hgjb1_second[0:50])
# coordinates1.extend(coords_hgjb1_second[51:100])
# coordinates1.extend(coords_hgjb1_second[101:150])
# coordinates1.extend(coords_hgjb1_second[151:200])

# edge_hgjb1_second = cq.Edge.makeSpline([cq.Vector(p) for p in coordinates1][0:199])





# edge_hgjb1_first = edge_hgjb1_first_1.fuse(edge_hgjb1_first_2, glue = True).clean().fix()
# edge_hgjb1_first = edge_hgjb1_first.fuse(edge_hgjb1_first_3, glue = True).clean().fix()
# edge_hgjb1_first = edge_hgjb1_first.fuse(edge_hgjb1_first_4, glue = True).clean().fix()

# show_object(edge_hgjb1_first)
# show_object(edge_hgjb1_second)



# edge_hgjb1_second = cq.Edge.makeSpline([cq.Vector(p) for p in coords_hgjb1_second][0:points]).close()
# show_object(edge_hgjb1_second)
# edge_hgjb2_first = cq.Edge.makeSpline([cq.Vector(p) for p in coords_hgjb2_first][0:points]).close()
# show_object(edge_hgjb2_first)
# edge_hgjb2_second = cq.Edge.makeSpline([cq.Vector(p) for p in coords_hgjb2_second][0:points]).close()
# show_object(edge_hgjb2_second)

# face_hgjb1_first = cq.Face.makeRuledSurface(edge_hgjb1_first_1,edge_hgjb1_first_3)
# show_object(face_hgjb1_first)


# face_hgjb1_first_1 = cq.Face.makeRuledSurface(edge_hgjb1_first_5,edge_hgjb1_first_6)
# show_object(face_hgjb1_first_1)

# rule1 = face_hgjb1_first.fuse(face_hgjb1_first_1, glue = True).clean().fix()


# loft_hgjb1 = cq.Solid.makeLoft([edge_hgjb1_first,edge_hgjb1_second],True)



# face_hgjb1_second = cq.Face.makeRuledSurface(edge_hgjb1_second_1,edge_hgjb1_second_3)
# show_object(face_hgjb1_second)



# face_hgjb1_first = cq.Face.makeNSidedSurface([edge_hgjb1_first],[],tol2d=0.0001,tol3d=0.0001,tolAng=0.0001,tolCurv=0.00001,maxDeg=2,maxSegments=2,nbIter=3)
# show_object(face_hgjb1_first)

# face_hgjb1_second = cq.Face.makeNSidedSurface([edge_hgjb1_second],[],tol2d=0.0001,tol3d=0.0001,tolAng=0.0001,tolCurv=0.00001,maxDeg=2,maxSegments=222,nbIter=3)
# show_object(face_hgjb1_second)



# face_hgjb1_second = cq.Face.makeNSidedSurface([edge_hgjb1_second],[])
# show_object(face_hgjb1_second)


# face_hgjb2_first = cq.Face.makeNSidedSurface([edge_hgjb2_first],[])
# face_hgjb2_second = cq.Face.makeNSidedSurface([edge_hgjb2_second ],[])

# loft_hgjb1 = cq.Solid.makeLoft([edge_hgjb1_first,edge_hgjb1_second],True)
# loft_hgjb2 = cq.Solid.makeLoft([edge_hgjb2_first, edge_hgjb2_second],True)

# shell_hgjb1 = cq.Shell.makeShell([face_hgjb1_first, face_hgjb1_second, loft_hgjb1]).fix()
# shell_hgjb2 = cq.Shell.makeShell([face_hgjb2_first, face_hgjb2_second, loft_hgjb2]).fix()

# #calculate turn angle in radians
# radang = (gap/CylRadOut1)#gap/CylRadOut1 #gap_spiral/(pi*D)# Betaprime #gap_spiral/pi*D #gap/CylRadOut1

# #convert to degrees
# rotang = radang*180/pi 

# solid_hgjb1 ={}
# solid_hgjb1[0] = cq.Solid.makeSolid(shell_hgjb1)
# solid_hgjb1[0] = solid_hgjb1[0].transformed((0,0,0),(0,-7.5,10))
# show_object(solid_hgjb1[0])

# solid_hgjb2 ={}
# solid_hgjb2 [0] = cq.Solid.makeSolid(shell_hgjb2)


# for i in range(0,28):
#     solid_hgjb1[i+1] = solid_hgjb1[i].transformed ((0 ,0 ,sepang))
#     #show_object(solid_hgjb1[i+1])
#     #Rotor = Rotor.cut(solid_hgjb1[i+1])
    
#     # solid_hgjb2[i+1] = solid_hgjb2[i].transformed ((0 ,sepang ,0))
#     # Rotor = Rotor.cut(solid_hgjb2[i+1])
    
#     # solid_hgjb1m[i+1] = solid_hgjb1m[i].transformed ((0 ,sepang ,0))
#     # Rotor = Rotor.cut(solid_hgjb1m[i+1])
    

#show_object(Rotor)

#cq.exporters.export(cylinder1, "HGJBonRotor_Element_23_08_19.step")






