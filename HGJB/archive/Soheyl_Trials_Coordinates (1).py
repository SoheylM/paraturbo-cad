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

from cadquery import Edge, Vector, Wire, Solid,Shell

import time
import sys

from cqmore import Workplane
from cqmore.polyhedron import gridSurface

from cqmore import Workplane
from cqmore.polyhedron import gridSurface

from cqmore import Workplane
from cqmore.matrix import translation, rotationX, rotationZ
from cqmore.polyhedron import sweep
from cqmore.matrix import mirror, translation

t0 = time.time()

# Load Element
Element = pickle.load(open("Z:/Code/paraturbo-cad/ELEMENT/Element_23_06_05_v2.pickle", "rb"))

# Load the Rotor for the grooving test
Rotor = cq.importers.importStep("Z:/Code/paraturbo-cad/STEP/Rotor_23_06_05_v2.step").translate((0,0,-50))

# Capture main dimensions
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

inner_1 = Element['parameters']['hgjb2']['x_first_curve'][0:50]
inner_2 = Element['parameters']['hgjb2']['y_first_curve'][0:50]
inner_3 = Element['parameters']['hgjb2']['z_first_curve'][0:50]
inner1 = list(zip(inner_1,inner_2,inner_3))

inner_11 = Element['parameters']['hgjb2']['x_first_curve'][100:150]
inner_21 = Element['parameters']['hgjb2']['y_first_curve'][100:150]
inner_31 = Element['parameters']['hgjb2']['z_first_curve'][100:150]
inner11 = list(zip(inner_11,inner_21,inner_31))


x_second = Element['parameters']['hgjb2']['x_second_curve'][0:50]
y_second = Element['parameters']['hgjb2']['y_second_curve'][0:50]
z_second = Element['parameters']['hgjb2']['z_second_curve'][0:50]
inner2 = list(zip(x_second,y_second,z_second))

x_second1 = Element['parameters']['hgjb2']['x_second_curve'][100:150]
y_second1 = Element['parameters']['hgjb2']['y_second_curve'][100:150]
z_second1 = Element['parameters']['hgjb2']['z_second_curve'][100:150]
inner22 = list(zip(x_second1,y_second1,z_second1))



inner11.reverse()

inner22.reverse()


inner = [inner1,inner11]



# surface = Workplane().polyhedron(*gridSurface(inner,0.5))

# surface1 = Workplane().splineApproxSurface(inner,0.1)


# # show_object(surfacee)

# shell = Shell.makeShell([surface1]).fix()

# rr = Solid.makeSolid(shell)
# show_object(rr)


# surface1 = Workplane().splineApproxSurface(inner,0.1)
# show_object(surface1)


# shell = cq.Shell.makeShell([surface1]).fix()


# solid = cq.Solid.makeSolid(surface1)

# Rotor = Rotor.cut(surface1)

# show_object(Rotor)



profile = {}
profile1 = {}

usage={}
usagemm={}
surface={}

usa={}
mir={}
mir1={}
surf={}

usa2={}
mir2={}
mir12={}
surf2={}


surfacemm={}

mirrored_pts={}
mirrored_pts1={}

trans={}
trans1={}
tra={}
tra1={}

angle_step = 360/28
profile[0] = inner1
profile1[0] = inner11
usage[0]=[inner1,inner11]
# surface[0] = Workplane().polyhedron(*gridSurface(inner,0.1))

# Rotor = Workplane().cylinder(500,7.42)

trans[0] = translation((0, 0, 1)).transformAll(profile[0])
trans1[0] = translation((0, 0, 1)).transformAll(profile1[0])

mirrored_pts[0] = mirror((0, 0, 1)).transformAll(trans[0])
mirrored_pts1[0] = mirror((0, 0, 1)).transformAll(trans1[0])
usagemm[0]=[mirrored_pts[0],mirrored_pts1[0]]

dist_hack = 60 #60
mir[0] = translation((0, 0, dist_hack)).transformAll(profile[0])
mir1[0] = translation((0, 0, dist_hack)).transformAll(profile1[0])
usa[0]=[mir[0],mir1[0]]

tra[0] = translation((0, 0, dist_hack)).transformAll(mirrored_pts[0])
tra1[0] = translation((0, 0, dist_hack)).transformAll(mirrored_pts1[0])
usa2[0]=[tra[0],tra1[0]]


mir2[0] = mirror((0, 0, -1)).transformAll(tra[0])
mir12[0] = mirror((0, 0, -1)).transformAll(tra1[0])
usa2[0]=[mir2[0],mir12[0]]

m = rotationZ(angle_step)


for i in range(0,28):
    profile[i+1]=m.transformAll(profile[i])
    profile1[i+1]=m.transformAll(profile1[i])
    usage[i+1]=[profile[i+1],profile1[i+1]]
    # surface[i+1] = Workplane().polyhedron(*gridSurface(usage[i],0.1))
    # 
    surface[i+1] = Workplane().splineApproxSurface(usage[i],-0.1,clean=True,combine=False)
    #show_object(surface[i+1])
    # Rotor = Rotor - surface[i+1]

    #trans[i+1] = translation((0, 0, 1)).transformAll(profile[i])
    #trans1[i+1] = translation((0, 0, 1)).transformAll(profile1[i])

    #mirrored_pts[i+1] = mirror((0, 0, 1)).transformAll(trans[i])
    #mirrored_pts1[i+1] = mirror((0, 0, 1)).transformAll(trans1[i])

    mirrored_pts[i+1]=m.transformAll(mirrored_pts[i])
    mirrored_pts1[i+1]=m.transformAll(mirrored_pts1[i])
    usagemm[i+1]=[mirrored_pts[i+1],mirrored_pts1[i+1]]
    # surfacemm[i+1] = Workplane().polyhedron(*gridSurface(usagemm[i],0.1))
    # 
    surfacemm[i+1] = Workplane().splineApproxSurface(usagemm[i],0.1,clean=True,combine=False)
    # show_object(surfacemm[i+1])
    # Rotor = Rotor - surfacemm[i+1]

    mir[i+1]=m.transformAll(mir[i])
    mir1[i+1]=m.transformAll(mir1[i])
    usa[i+1]=[mir[i+1],mir1[i+1]]
    # surf[i+1] = Workplane().polyhedron(*gridSurface(usa[i],0.1))
    # 
    # show_object(surf[i+1])
    surf[i+1] = Workplane().splineApproxSurface(usa[i],-0.1,clean=True,combine=False)
    # Rotor = Rotor - surf[i+1]

    tra[i+1]=m.transformAll(tra[i])
    tra1[i+1]=m.transformAll(tra1[i])
    usa2[i+1]=[tra[i+1],tra1[i+1]]
    # surf2[i+1] = Workplane().polyhedron(*gridSurface(usa2[i],0.1))
    # Rotor = Rotor - surf2[i+1]
    # show_object(surf2[i+1])
    surf2[i+1] = Workplane().splineApproxSurface(usa2[i],0.1,clean=True,combine=False)
    # Rotor = Rotor - surf2[i+1]



# for i in range(0,3):

#     mir2[i+1]=m.transformAll(mir2[i])
#     mir12[i+1]=m.transformAll(mir12[i])
#     usa2[i+1]=[mir2[i+1],mir12[i+1]]
#     surf2[i+1] = Workplane().polyhedron(*gridSurface(usa2[i],0.1))
#     show_object(surf2[i+1])

# assembly = cq.Assembly()
# for i in range(28):
#     assembly.add(surface[i+1])
#     assembly.add(surfacemm[i+1])
#     assembly.add(surf[i+1])
#     assembly.add(surf2[i+1])

#assembly.save(self.cwf  + '/STEP/Compressor.step')

texts = Workplane()
for i in range(0,28): # to include the initial groove
    texts.add(surface[i+1])
    texts.add(surfacemm[i+1])
    texts.add(surf[i+1])
    texts.add(surf2[i+1])

# cq.exporters.export(texts,"grooves.step")
# grooves = cq.importers.importStep("C:/Users/user/Documents/Code/CadQuery/Compressor JB/paraturbo-cad/grooves.step")

# show_object(grooves)
    
# show_object(texts)

# show_object(assembly)
   
    
# Rotor = Rotor.cut(texts)

# show_object(texts)

Rotordone = Rotor - texts

# Rotordone = Workplane().cylinder(500,7.42) - texts
print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')

# show_object(Rotordone)

# print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')
# cq.exporters.export(Rotordone,"Rotordone.stl", tolerance = 0.1, angularTolerance = 0.5)
#cq.exporters.export(Rotordone,"Rotordone.stl", tolerance = 0.1, angularTolerance = 0.5)

#print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')
cq.exporters.export(Rotordone,"Rotordone.step", opt={"write_pcurves": False, "precision_mode": 1})
print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')




'''



# texts = texts.intersect(Rotor)
# show_object(Rotordone)
# print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')




# polyline = (Workplane()
            #     .polylineJoin(
            #         inner1, 
            #         Workplane().box(.1, .1, .1)
            #     )
            # )
# show_object(polyline)

# r = Workplane().polyhedron(*sweep(inner,closeIdx = 2))

# show_object(r)




# edge = cq.Edge.makeSpline([cq.Vector(p) for p in inner1])
# edge1 = cq.Edge.makeSpline([cq.Vector(p) for p in inner11])

# show_object(edge)
# show_object(edge1)

#show_object(Rotor)

# surface1 = Workplane().splineApproxSurface(inner,0.1)
# show_object(surface1)


#shell = cq.Shell.makeShell([surface1]).fix()


#solid = cq.Solid.makeSolid(surface1)

# Rotor = Rotor.cut(surface1)

# show_object(Rotor)


#solid ={}
# solid[0] = surface1


# for i in range(0,28):
#     solid[i+1] = solid[i].transformed((0,0,sepang))
#     show_object(solid[i+1])
#     Rotor = Rotor.cut(solid[i+1])
    
#     # solid_hgjb2[i+1] = solid_hgjb2[i].transformed ((0 ,sepang ,0))
#     # Rotor = Rotor.cut(solid_hgjb2[i+1])
    
#     # solid_hgjb1m[i+1] = solid_hgjb1m[i].transformed ((0 ,sepang ,0))
#     # Rotor = Rotor.cut(solid_hgjb1m[i+1])
    

# show_object(Rotor)

# profile = inner1
# profile1 = inner11

'''
# inner = [inner1,inner11]

# profile = {}
# profile1 = {}
# usage={}
# surface={}

# angle_step = 360/28
# profile[0] = inner1
# profile1[0] = inner11
# usage[0]=[inner1,inner11]
# surface[0] = Workplane().splineApproxSurface(usage[0],0.1)

# for i in range(0,2):
#     m = rotationZ(angle_step)
#     profile[i+1]=m.transformAll(profile[i])
#     profile1[i+1]=m.transformAll(profile1[i])
#     usage[i+1]=[profile[i+1],profile1[i+1]]
#     surface[i+1] = Workplane().splineApproxSurface(usage[i],0.1)
#     show_object(surface[i+1])
#     Rotor = Rotor.cut(surface[i+1])

# show_object(Rotor)
'''


# profile2 = {}
# profile12 = {}
# usage2={}
# surface2={}

# angle_step = 360/28
# profile2[0] = inner2
# profile12[0] = inner22
# usage2[0]=[inner2,inner22]
# surface2[0] = Workplane().splineApproxSurface(usage2[0],0.01)

# for i in range(0,28):
#     m = rotationZ(angle_step)
#     profile2[i+1]=m.transformAll(profile2[i])
#     profile12[i+1]=m.transformAll(profile12[i])
#     usage2[i+1]=[profile2[i+1],profile12[i+1]]
#     surface2[i+1] = Workplane().splineApproxSurface(usage2[i],0.01)
#     show_object(surface2[i+1])





















# inner_111 = Element1['parameters']['hgjb1']['x_first_surface'][101:150]
# inner_211 = Element1['parameters']['hgjb1']['y_first_surface'][101:150]
# inner_311 = Element1['parameters']['hgjb1']['z_first_surface'][101:150]
# inner111 = list(zip(inner_111,inner_211,inner_311))



# inner = [coords_hgjb1_first]
'''
# hgjb1_x_second = Element['parameters']['hgjb1']['x_second_curve']
# hgjb1_y_second = Element['parameters']['hgjb1']['y_second_curve']
# hgjb1_z_second = Element['parameters']['hgjb1']['z_second_curve']
# coords_hgjb1_second = list(zip(hgjb1_x_second, hgjb1_y_second, hgjb1_z_second))

# hgjb2_x_first = Element['parameters']['hgjb2']['x_first_curve']
# hgjb2_y_first = Element['parameters']['hgjb2']['y_first_curve']
# hgjb2_z_first = Element['parameters']['hgjb2']['z_first_curve']
# coords_hgjb2_first = list(zip(hgjb2_x_first, hgjb2_y_first, hgjb2_z_first))

# hgjb2_x_second = Element['parameters']['hgjb2']['x_second_curve']
# hgjb2_y_second = Element['parameters']['hgjb2']['y_second_curve']
# hgjb2_z_second = Element['parameters']['hgjb2']['z_second_curve']
# coords_hgjb2_second = list(zip(hgjb2_x_second, hgjb2_y_second, hgjb2_z_second))



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

# cq.exporters.export(Rotor,"JBtrial.step")




#################
###################
#####################


# from cqmore import Workplane

# def paraboloid(x, y):
#     return (x, y, ((y ** 2) - (x ** 2)) / 4)

# min_value = -30
# max_value = 30
# step = 5
# thickness = 0.5

# points = [[
#         paraboloid(x / 10, y / 10) 
#     for y in range(min_value, max_value + step, step)
# ] for x in range(min_value, max_value + step, step)]

# surface33 = Workplane().splineApproxSurface(points, thickness)

# show_object(surface33)





#################
###################
#####################