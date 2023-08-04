import pickle
import numpy as np
import cadquery as cq
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions

#creating the rotor

#open file
file = open('best_solution_Element.pickle', 'rb')

# dump info to that file
Element = pickle.load(file)

#close file
file.close()

#print('pos_hgjb1',Element.pos_hgjb1)

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

sys_pos = Element['sys_pos']
pos_hgjb1 = sys_pos['pos_hgjb1']
pos_hgjb2 =sys_pos['pos_hgjb2']
print(pos_hgjb1)
print(pos_hgjb2)

#NOTE: circle is defined by radius
#NOTE: hole is defined by diameter

WP = cq.Workplane('XY')

#part=WP.cylinder(10,2, centered =(True,True,False)) centers the cylinder on x and y axis but not on z axis

#initialize part
#start outside in: draw outermost cylinder, then make a hole
#this avoids drawing a cylinder, drawing another cylinder around it and putting a hole through both

#test presence of outermost cylinder, if found, draw and test presence of middle cylinder and inner cylinder, draw if present
if DA3[0] != DI3[0]:
    part3 = WP.cylinder(Laenge[0],DA3[0]/2, direct=(0,0,1), angle = 360, centered=(True,True,False),combine=False).hole(DI3[0],Laenge[0])
    show_object(part3)
    if DA2[0] != DI2[0]:
        part2=WP.cylinder(Laenge[0],DA2[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI2[0],Laenge[0])    
        show_object(part2)
        if DA1[0] != DI1[0]:
                part1=WP.cylinder(Laenge[0],DA1[0]/2,direct=(0,0,1), angle = 360, centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
                show_object(part1)
    elif DA1[0] != DI1[0]:
            part1=WP.cylinder(Laenge[0],DA1[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
            show_object(part1)
#if no outermost cylinder, test presence of middle cyinder, draw if present, test for innermost cylinder and draw if present
elif DA2[0] != DI2[0]:
    part2=WP.cylinder(Laenge[0],DA2[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI2[0],Laenge[0])
    show_object(part2)
    if DA1[0] != DI1[0]:
        part1=WP.cylinder(Laenge[0],DA1[0]/2,direct=(0,0,1), angle = 360, centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
        show_object(part1)
#if no outermost and middle cylinder, test for inner cylinder and draw if present
elif DA1[0] != DI1[0]:
    part1 =WP.cylinder(Laenge[0],DA1[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
    show_object(part1)
#if no cylinders present, print informational message
else:
    print("First section empty")
#part=WP.cylinder(Laenge[0],DA1[0]/2, centered=(True,True,False)).hole(DI1[0],Laenge[0])

#build successive cylinder stacks
#for outermost cylinder present, move to uppermost z face and add next cylinder
#check that the right face is being used at each stacking
for i in range(21):
    if DA3[i] != DI3[i]:
        WP = part3.faces('>Z')
    elif DA2[i] != DI2[i]:
        WP = part2.faces('>Z')
    elif DA1[i] != DI1[i]:
        WP = part1.faces('>Z')

    if DA3[i+1] != DI3[i+1]:
        part3=WP.cylinder(Laenge[i+1],DA3[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI3[i+1],Laenge[i+1])
        show_object(part3)
        if DA2[i+1] != DI2[i+1]:
            part2=WP.cylinder(Laenge[i+1],DA2[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI2[i+1],Laenge[i+1])
            show_object(part2)
        elif DA1[i+1] != DI1[i+1]:
            part1=WP.cylinder(Laenge[i+1],DA1[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI1[i+1],Laenge[i+1])
            show_object(part1)
    elif DA2[i+1] != DI2[i+1]:
        #print('DA2[i+1] != DI2[i+1]',i+1)
        part2=WP.cylinder(Laenge[i+1],DA2[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI2[i+1],Laenge[i+1])
        show_object(part2)
        if DA1[i+1] != DI1[i+1]:
            #print('DA1[i+1] != DI1[i+1]',i+1)
            part1=WP.cylinder(Laenge[i+1],DA1[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI1[i+1],Laenge[i+1])
            show_object(part1)
    elif DA1[i+1] != DI1[i+1]:
        part1=WP.cylinder(Laenge[i+1],DA1[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI1[i+1],Laenge[i+1])
        show_object(part1)



#creating the grooves

#number of units present
units=28

#calculate angle of separation based on how many units present
sepang=360/units

#Length of entire bearing along cylinder axis
Len = Laenge[pos_hgjb1]

# #Length of middle strip
L_land = 3.08 #mm

# #Distance to center of bearing along rotor axis
for d in range(pos_hgjb1):
    
# DistCenter = Len/2 #+sum of all other cylinder lengths

# #Outside angle between parallelogram vertical and diagonal
# Beta = 135

# #Groove height (parallelogram vertical)
# a = 1

# #length between parallelogram verticals
# LenBetwVert = Len/2 - Lland/2

# #Angle between parallelogram diagonal and reference horizontal
# Betaprime = Beta - 90

# #gap between horizontal leaving from parallelogram lowest 
# #corner and parallelogram vertical

# gap = LenBetwVert*tan(Betaprime)


# #Length of cylinder section
# CylLen = 100

# #Outer radius of the cylinder
# CylRadOut = 50

# #cylinder to be projected onto (100 is the length)
# cylinder = cq.Solid.makeCylinder(
#      CylRadOut, CylLen, pnt=cq.Vector(0, -50, 0), dir=cq.Vector(0, 1, 0)
# )

# #cylinder = cylinder.transformed(rotate=(0,60,0))

# removalcylinder = cq.Solid.makeCylinder(
#       CylRadOut, CylLen, pnt=cq.Vector(0, -50, -30), dir=cq.Vector(0, 1, 0)
# )

# #direction of projection 
# projection_direction = cq.Vector(0, 0, 1)

# yp1 = DistCenter 
# xp1 = LenBetwVert*tan(Betaprime)/2 
# yp2 = -LenBetwVert*tan(Betaprime)
# xp2 = -LenBetwVert

# #draw first parallelogram as a 3D shape in YX plane
# parallelogram1 = (
#       cq.Workplane("YX", origin=((gap+a)/2, yp1-1, -2*CylRadOut))
#       #when viewed / \, y to the right and x up, z into screen
#       #origin at upper outside corner
#       #points below in standard x and y coordinates
#       .lineTo(-LenBetwVert,-gap) #upper inside corner 
#       .lineTo(-LenBetwVert,-gap-a) #lower inside corner
#       .lineTo(0,-a) #lower outside corner
#       .close()
#       .extrude(1)
#       .faces("<Z")
#       .val()
#   )

# parallelogram2 = (
#       cq.Workplane("YX", origin=((gap+a)/2, -yp1+1, -2*CylRadOut))
#       #when viewed / \, y to the right and x up, z into screen
#       #origin at upper outside corner
#       #points below in standard x and y coordinates
#       .lineTo(LenBetwVert,-gap) #upper inside corner 
#       .lineTo(LenBetwVert,-gap-a) #lower inside corner
#       .lineTo(0,-a) #lower outside corner
#       .close()
#       .extrude(1)
#       .faces("<Z")
#       .val()
#   )

# # #height between outer top and inner top
# # h=60
# # #vertical height
# # v=1

# # parallelogram1 = (
# #       cq.Workplane("YX", origin=((h+v)/2, 40, -2*CylRadOut))
# #       #when viewed / \, y to the right and x up, z into screen
# #       #origin at upper outside corner
# #       .lineTo(-40,-h) #upper inside corner 
# #       .lineTo(-40,-h-v) #lower inside corner
# #       .lineTo(0,-v) #lower outside corner
# #       .close()
# #       .extrude(1)
# #       .faces("<Z")
# #       .val()
# #   )


# #project first parallelogram onto cylinder
# parallelogram1_projected = parallelogram1.projectToShape(cylinder, projection_direction)

# #turn first parallelogram into 3D shape on cylinder surface
# parallelogram1_solids = cq.Compound.makeCompound(
#       [f.thicken(2) for f in parallelogram1_projected]
#   )
# parallelogram1_solids = parallelogram1_solids.cut(removalcylinder)

# # #project second parallelogram onto cylinder
# # parallelogram2_projected = parallelogram2.projectToShape(cylinder, projection_direction)

# # #turn second parallelogram into 3D shape on cylinder surface
# # parallelogram2_solids = cq.Compound.makeCompound(
# #      [f.thicken(2) for f in parallelogram2_projected]
# #  )

# # parallelogram2_solids = parallelogram2_solids.cut(removalcylinder)

# # #removes the projected parallelogram 3D objects from the cylinder

# # for i in range(units):
# #     cylinder = cylinder.cut(parallelogram1_solids)
# #     #cylinder = cylinder.cut(parallelogram2_solids)
# #     cylinder = cylinder.transformed(rotate=(0,sepang,0))

# # cylinder=cylinder.transformed(rotate=(0,0,0))

# #show the various objects
# if "show_object" in locals():
#      show_object(cylinder, name="cylinder_solid", options={"alpha": 0.8})
#      #show_object(removalcylinder)
#      show_object(parallelogram1, name="parallelogram1")
#      show_object(parallelogram2, name="parallelogram2")
     
#      #these are the two objects we are interested in
#      #each one creates a curved parallelogram outside and
#      #inside the cylinder; need to show both and remove excess
#      #show_object(parallelogram1_solids, name="parallelogram1_solids")
#      #show_object(parallelogram2_solids, name="parallelogram2_solids")



