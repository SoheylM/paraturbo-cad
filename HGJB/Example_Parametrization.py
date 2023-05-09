from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cadquery as cq
import cq_warehouse.extensions

#number of units present
units=20

#calculate angle of separation based on how many units present
sepang=360/units

#Length of entire bearing along cylinder axis
Len = 100

#Length of middle strip
Lland = 3

#Distance to center of bearing along rotor axis
DistCenter = Len/2 #+sum of all other cylinder lengths

#Outside angle between parallelogram vertical and diagonal
Beta = 135

#Groove height (parallelogram vertical)
a = 1

#length between parallelogram verticals
LenBetwVert = Len/2 - Lland/2

#Angle between parallelogram diagonal and reference horizontal
Betaprime = Beta - 90

#gap between horizontal leaving from parallelogram lowest 
#corner and parallelogram vertical

gap = LenBetwVert*tan(Betaprime)


#Length of cylinder section
CylLen = 100

#Outer radius of the cylinder
CylRadOut = 50

#cylinder to be projected onto (100 is the length)
cylinder = cq.Solid.makeCylinder(
     CylRadOut, CylLen, pnt=cq.Vector(0, -50, 0), dir=cq.Vector(0, 1, 0)
)

#cylinder = cylinder.transformed(rotate=(0,60,0))

removalcylinder = cq.Solid.makeCylinder(
      CylRadOut, CylLen, pnt=cq.Vector(0, -50, -50), dir=cq.Vector(0, 1, 0)
)

#direction of projection 
projection_direction = cq.Vector(0, 0, 1)

yp1 = DistCenter 
xp1 = LenBetwVert*tan(Betaprime)/2 
yp2 = -LenBetwVert*tan(Betaprime)
xp2 = -LenBetwVert

#draw first parallelogram as a 3D shape in YX plane
parallelogram1 = (
     cq.Workplane("YX", origin=((gap+a)/2, yp1, -2*CylRadOut))
     #when viewed / \, y to the right and x up, z into screen
     #origin at upper outside corner
     #points below in standard x and y coordinates
     .lineTo(-LenBetwVert,-gap) #upper inside corner 
     .lineTo(-LenBetwVert,-gap-a) #lower inside corner
     .lineTo(0,-a) #lower outside corner
     .close()
     .extrude(1)
     .faces("<Z")
     .val()
 )

# parallelogram1 = (
#      cq.Workplane("YX", origin=(10, 30, -2*CylRadOut))
#      #when viewed / \, y to the right and x up, z into screen
#      #origin at upper outside corner
#      .lineTo(-30,10) #upper inside corner 
#      .lineTo(-30,0) #lower inside corner
#      .lineTo(0,-10) #lower outside corner
#      .close()
#      .extrude(1)
#      .faces("<Z")
#      .val()
#  )

#draw second parallelogram as a 3D shape in YX plane
# parallelogram2 = (
#      cq.Workplane("YX", origin=(0, -35, -2*CylRadOut))
#      #when viewed / \, y to the left and x down, z into screen
#      #origin at upper outside corner
#      .lineTo(30,10) #upper inside corner
#      .lineTo(30,0) #lower inside corner
#      .lineTo(0,-10) #lower outside corner
#      .close()
#      .extrude(1)
#      .faces("<Z")
#      .val()
#  )
# #project first parallelogram onto cylinder
# parallelogram1_projected = parallelogram1.projectToShape(cylinder, projection_direction)

# #turn first parallelogram into 3D shape on cylinder surface
# parallelogram1_solids = cq.Compound.makeCompound(
#      [f.thicken(2) for f in parallelogram1_projected]
#  )
# parallelogram1_solids = parallelogram1_solids.cut(removalcylinder)

# #project second parallelogram onto cylinder
# parallelogram2_projected = parallelogram2.projectToShape(cylinder, projection_direction)

# #turn second parallelogram into 3D shape on cylinder surface
# parallelogram2_solids = cq.Compound.makeCompound(
#      [f.thicken(2) for f in parallelogram2_projected]
#  )

# parallelogram2_solids = parallelogram2_solids.cut(removalcylinder)

# #removes the projected parallelogram 3D objects from the cylinder

# for i in range(units):
#     cylinder = cylinder.cut(parallelogram1_solids)
#     cylinder = cylinder.cut(parallelogram2_solids)
#     cylinder = cylinder.transformed(rotate=(0,sepang,0))

# cylinder=cylinder.transformed(rotate=(0,0,0))

#show the various objects
if "show_object" in locals():
     show_object(cylinder, name="cylinder_solid") #, options={"alpha": 0.8}
     #show_object(removalcylinder)
     show_object(parallelogram1, name="parallelogram1")
     #show_object(parallelogram2, name="parallelogram2")
     
     #these are the two objects we are interested in
     #each one creates a curved parallelogram outside and
     #inside the cylinder; need to show both and remove excess
     #show_object(parallelogram1_solids, name="parallelogram1_solids")
     #show_object(parallelogram2_solids, name="parallelogram2_solids")
