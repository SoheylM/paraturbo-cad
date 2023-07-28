#adapted from May 11 version of takingrotor.py

import pickle
import numpy as np
import cadquery as cq
from math import pi, tan
#from cq_warehouse import *
import time
from cqmore import Workplane
from cqmore.matrix import translation, rotationZ
from cqmore.matrix import mirror, translation

t0 = time.time()

# Load Element
Element = pickle.load(open("Z:/Code/paraturbo-cad/ELEMENT/Element_23_06_05_v2.pickle", "rb"))


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
print('L_hgjb', L)
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
print('DistCenter1', DistCenter1)

#find distance to center of second HGJB
DistCenter2=dist2+L/2
print('DistCenter2', DistCenter2)


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



######### START OF IMPLEMENTATION ##########
# Length to translate hgjb1
dist_hack0 = -(DistCenter1+0-L_land/2) #-50
# Length to increase/decrease L_land
dist_hack1 = 0 #1
# Length to translate hgjb2
dist_hack2 = DistCenter2 -(DistCenter1+0) #60


# Load the Rotor for the grooving test
Rotor = cq.importers.importStep("Z:/Code/paraturbo-cad/STEP/Rotor_23_06_05_v2.step").translate((0,0,dist_hack0))


# Extract grooves coordinates from Element
inner_1 = Element['parameters']['hgjb2']['x_first_curve'][0:50]
inner_2 = Element['parameters']['hgjb2']['y_first_curve'][0:50]
inner_3 = Element['parameters']['hgjb2']['z_first_curve'][0:50]
inner1 = list(zip(inner_1,inner_2,inner_3))

inner_11 = Element['parameters']['hgjb2']['x_first_curve'][100:150]
inner_21 = Element['parameters']['hgjb2']['y_first_curve'][100:150]
inner_31 = Element['parameters']['hgjb2']['z_first_curve'][100:150]
inner11 = list(zip(inner_11,inner_21,inner_31))
inner11.reverse()

# Initialize dictionaries
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

angle_step  = 360/N_HG
profile[0]  = inner1
profile1[0] = inner11
usage[0]    = [inner1,inner11]

#dist_hack1 = 1 #1
trans[0]  = translation((0, 0, dist_hack1)).transformAll(profile[0])
trans1[0] = translation((0, 0, dist_hack1)).transformAll(profile1[0])

mirrored_pts[0]  = mirror((0, 0, 1)).transformAll(trans[0])
mirrored_pts1[0] = mirror((0, 0, 1)).transformAll(trans1[0])
usagemm[0]       = [mirrored_pts[0],mirrored_pts1[0]]

#dist_hack2 = 60 #60
mir[0]    = translation((0, 0, dist_hack2)).transformAll(profile[0])
mir1[0]   = translation((0, 0, dist_hack2)).transformAll(profile1[0])
usa[0]    = [mir[0],mir1[0]]

tra[0]    = translation((0, 0, dist_hack2)).transformAll(mirrored_pts[0])
tra1[0]   = translation((0, 0, dist_hack2)).transformAll(mirrored_pts1[0])
usa2[0]   = [tra[0],tra1[0]]

m = rotationZ(angle_step)

#h_gr = 1 #overwrite

for i in range(0,28):
    profile[i+1]  = m.transformAll(profile[i])
    profile1[i+1] = m.transformAll(profile1[i])
    usage[i+1]    = [profile[i+1],profile1[i+1]]
    surface[i+1]  = Workplane().splineApproxSurface(usage[i],-0.1,clean=True,combine=False)
    if i == -1:
        print('usage 0', usage[i])

    mirrored_pts[i+1]  = m.transformAll(mirrored_pts[i])
    mirrored_pts1[i+1] = m.transformAll(mirrored_pts1[i])
    usagemm[i+1]       = [mirrored_pts[i+1],mirrored_pts1[i+1]]
    surfacemm[i+1]     = Workplane().splineApproxSurface(usagemm[i],0.1,clean=True,combine=False)
    if i == -1:
        print('usagemm 0', usagemm[i])

    mir[i+1]=m.transformAll(mir[i])
    mir1[i+1]=m.transformAll(mir1[i])
    usa[i+1]=[mir[i+1],mir1[i+1]]
    surf[i+1] = Workplane().splineApproxSurface(usa[i],-0.1,clean=True,combine=False)
    if i == -1:
        print('usa 0', usa[i])

    tra[i+1]=m.transformAll(tra[i])
    tra1[i+1]=m.transformAll(tra1[i])
    usa2[i+1]=[tra[i+1],tra1[i+1]]
    surf2[i+1] = Workplane().splineApproxSurface(usa2[i],0.1,clean=True,combine=False)
    if i == -1:
        print('usa2 0', usa2[i])


texts = Workplane()
# each 'surface' is one of the half bearing, there are 4 in total
for i in range(0,28): 
    texts.add(surface[i+1])
    texts.add(surfacemm[i+1])
    texts.add(surf[i+1])
    texts.add(surf2[i+1])

Rotordone = Rotor - texts

print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')

cq.exporters.export(Rotordone,"Rotordone.step", opt={"write_pcurves": False, "precision_mode": 1})
print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')
