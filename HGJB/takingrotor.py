import pickle
import numpy as np
import cadquery as cq
import os
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions

#open file
file = open('Element_23_05_03.pickle', 'rb')

# dump info to that file
Element = pickle.load(file)

#close file
file.close()

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

cwf = os.getcwd().replace("\\", "/")

# Use this in CQ-Editor
Rotor = cq.importers.importStep(cwf + '/Rotor.stp')

# Use this in VS Code
# Rotor = cq.importers.importStep(cwf + '/HGJB/Rotor.stp')


Rotor = Rotor.rotate((0,0,0),(1,0,0),270)

#Begin code adapted from Christophe's Matlab
#Parameters 
N_HG = 28 #number of grooves generally between 26 - 30
alpha_HG = 0.68 #given
beta_HG = -135 #given 
beta_HG = beta_HG*pi/180
gamma_HG = 0.89 #given
h_gr = 16 #groove depth given in micrometers
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

#find distance to first HGJB
dist1 = 0
for d in range(pos_hgjb1):
    dist1=dist1+Laenge[d]
#find distance to center of first HGJB
DistCenter1=dist1+L/2

#define separation angle between grooves
sepang=360/N_HG

#length between parallelogram verticals
LenBetwVert = L/2 - L_land/2

#Angle between parallelogram diagonal and reference horizontal
Betaprime = abs(beta_HG) - 90

#gap between horizontal leaving from parallelogram lowest 
#corner and parallelogram vertical

gap = LenBetwVert*tan(Betaprime)

#create removal cylinder
CylLen = Laenge[pos_hgjb1]
CylRadOut= DA3[pos_hgjb1]/2

#cylinder to be projected onto
cylinder = cq.Solid.makeCylinder(
     CylRadOut, CylLen+2, pnt=cq.Vector(0, DistCenter1+L/2+1, 0), dir=cq.Vector(0, -1, 0)
)

removalcylinder = cq.Solid.makeCylinder(
      CylRadOut*2, CylLen, pnt=cq.Vector(0, DistCenter1+L/2, -DA3[pos_hgjb1]*0.8), dir=cq.Vector(0, -1, 0)
)

#direction of projection 
projection_direction = cq.Vector(0, 0, 1)

#global coordinates of origin of parallelogram1
yp1 = DistCenter1+L/2
xp1 = LenBetwVert*tan(Betaprime)/2 
yp2 = -LenBetwVert*tan(Betaprime)
xp2 = -LenBetwVert

#draw first parallelogram as a 3D shape in YX plane
parallelogram1 = (
      cq.Workplane("YX", origin=((gap+a_HG)/2, yp1, -2*CylRadOut))
      #when viewed / \, y to the right and x up, z into screen
      #origin at upper outside corner
      #points below in standard x and y coordinates
      .lineTo(-LenBetwVert,-gap) #upper inside corner 
      .lineTo(-LenBetwVert,-gap-a_HG) #lower inside corner
      .lineTo(0,-a_HG) #lower outside corner
      .close()
      .extrude(1)
      .faces("<Z")
      .val()
  )



#project first parallelogram onto cylinder
parallelogram1_projected = parallelogram1.projectToShape(cylinder, projection_direction)

#turn first parallelogram into 3D shape on cylinder surface
parallelogram1_solids = cq.Compound.makeCompound(
      [f.thicken(h_gr/1000) for f in parallelogram1_projected]
  )
parallelogram1_solids = parallelogram1_solids.cut(removalcylinder)

for i in range(N_HG):
    cylinder = cylinder.cut(parallelogram1_solids)
    #cylinder = cylinder.cut(parallelogram2_solids)
    cylinder = cylinder.transformed(rotate=(0,sepang,0))


#show_object(Rotor)
#show_object(removalcylinder)
show_object(cylinder)
show_object(parallelogram1, name="parallelogram1")
show_object(parallelogram1_solids, name="parallelogram1_solids")