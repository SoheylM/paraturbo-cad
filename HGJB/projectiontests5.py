#adapted from May 11 version of takingrotor.py

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

#create removal cylinder1
CylLen1 = Laenge[pos_hgjb1]
CylRadOut1= DA3[pos_hgjb1]/2


#first cylinder to be projected onto
cylinder1 = cq.Solid.makeCylinder(
     CylRadOut1, 2*CylLen1+2, pnt=cq.Vector(0, DistCenter1+2*L+1, 0), dir=cq.Vector(0, -1, 0)
)
#first removal cylinder
removalcylinder1 = cq.Solid.makeCylinder(
      CylRadOut1*2, CylLen1, pnt=cq.Vector(0, DistCenter1+2*L, -DA3[pos_hgjb1]*0.8), dir=cq.Vector(0, -1, 0)
)


#direction of projection 
projection_direction = cq.Vector(0, 0, -1)



#calculate turn angle in radians
radang = (gap/CylRadOut1)#gap/CylRadOut1 #gap_spiral/(pi*D)# Betaprime #gap_spiral/pi*D #gap/CylRadOut1

#convert to degrees
rotang = radang*180/pi #+0.273

#rotang = 24.5

parallelograms = []
parallelograms_projected = []
curve = cq.Assembly()


for i in range(n_parall):
    print('i=', i)
    parallelogram = (
          cq.Workplane("YX", origin=((gap+a_HG)/2, DistCenter1+L/2+i*LenBetwVert, 0))
          #when viewed / \, y to the right and x up, z into screen
          #origin at upper outside corner
          #points below in standard x and y coordinates
          .lineTo(-(LenBetwVert)*(1+eps_perc),-gap) #upper inside corner 
          .lineTo(-(LenBetwVert)*(1+eps_perc),-gap-a_HG) #lower inside corner
          .lineTo(0,-a_HG) #lower outside corner
          .close()
          .extrude(1)
          .faces("<Z")
          .val()
      )
    parallelograms.append(parallelogram)
    
for i in range(n_parall):
    if i == 0:
        cylinder1 = cylinder1.rotate((0,0,0),(0,1,0),rotang)
    parallelogram_projected = parallelograms[i].projectToShape(cylinder1, projection_direction)
    parallelograms_projected.append(parallelogram_projected)
        
for i in range(n_parall):        
    parallelogram_solids = cq.Compound.makeCompound(
        [f.thicken(1, cq.Vector(0, 0, 1)) for f in parallelograms_projected[i]]
        )

    
    #cuts rotates the cylinder all the way around and cuts
    # for j in range(N_HG):
    #     cylinder1 = cylinder1.cut(parallelogram_solids)
    #     cylinder1=cylinder1.rotate((0,0,0),(0,1,0),-1*sepang)
    
    
    cylinder1 = cylinder1.cut(parallelogram_solids)
    
    #rotates, cuts a second groove, rotates back
    # cylinder1=cylinder1.rotate((0,0,0),(0,1,0),-1*sepang)
    # cylinder1 = cylinder1.cut(parallelogram_solids)
    # cylinder1=cylinder1.rotate((0,0,0),(0,1,0),1*sepang)
    
    if i == n_parall-1:
        pass
    else:
        cylinder1 = cylinder1.rotate((0,0,0),(0,1,0),rotang)
        
# cylinder1=cylinder1.rotate((0,0,0),(0,1,0),-rotang*(n_parall-1))
# cylinder1=cylinder1.rotate((0,0,0),(0,1,0),-1*sepang)






#show_object(removalcylinder1)
show_object(cylinder1) #, options={"alpha": 0.8}
# show_object(parallelograms_projected[0])
show_object(parallelogram_solids)
