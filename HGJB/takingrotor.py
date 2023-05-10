import cadquery as cq
import os
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions

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
L = 28 #on drawing [mm]
L_land=L-(gamma_HG*L) #Value for CAD
L=L+0.8 #oversized length for safety generally between 0.6 - 1
Spiral_step = pi*D*tan(beta_HG)
Spiral_height = L/2
a_HG = (pi*D*alpha_HG)/N_HG #mm
a_HG_plus_b_HG = a_HG/alpha_HG #mm
h_rr_tot = h_rr*2 #diametral clearance given in micrometers

#end code adapted from Christophe's Matlab

show_object(Rotor)