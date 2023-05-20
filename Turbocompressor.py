import cadquery as cq
import numpy as np
import time
import os

from SGTB import SGTB
cwf = os.getcwd().replace("\\", "/")

import sys
sys.path.append(cwf  + '/SGTB')
sys.path.append(cwf  + '/ROT')
sys.path.append(cwf  + '/COMP')
sys.path.append(cwf  + '/HGJB')

from SGTB import *
from Rotor import *
from Impeller import *
from Helper import *

t0 = time.time()

'''
Example data for manual construction
'''

Length = [4.0, 3.0, 4.7, 6.0, 2.0, 8.011738903547421, 2.0049720014074204, 2.0079282926578075, 11.610315613390416, 33.74846795856666, 11.610315613390416, 2.005566980991642, 5.026995483891631, 16.0, 3.0]
DI1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 13.0, 0.0, 0.0]
DI2 = [11.0, 5.0, 5.0, 5.0, 13.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 13.0, 13.0, 13.0]
DI3 = [11.0, 5.0, 5.0, 5.0, 13.0, 13.0, 27.5, 13.0, 14.0, 14.0, 14.0, 13.899999999999999, 13.899999999999999, 13.899999999999999, 13.899999999999999]
DO1 = [11.0, 5.0, 5.0, 5.0, 13.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 13.0, 13.0, 13.0]
DO2 = [11.0, 5.0, 5.0, 5.0, 13.0, 13.0, 27.5, 13.0, 14.0, 14.0, 14.0, 13.899999999999999, 13.899999999999999, 13.899999999999999, 13.899999999999999]
DO3 = [11.0, 11.0, 15.6, 25.0, 13.0, 13.0, 27.5, 13.0, 14.0, 14.0, 14.0, 13.899999999999999, 13.899999999999999, 13.899999999999999, 13.899999999999999]
types1 = ['COMP1','COMP1','COMP1','COMP1','PLUG','PLUG','ROT','ROT','ROT','ROT','ROT','ROT','ROT','MAG','MAG']
types2 = ['COMP1','COMP1','COMP1','COMP1','PLUG','ROT','ROT','ROT','ROT','ROT','ROT','ROT','ROT','ROT','ROT']
types3 = ['COMP1','COMP1','COMP1','COMP1','PLUG','ROT','ROT','ROT','ROT','ROT','ROT','ROT','ROT','ROT','ROT']
pos = 6
alpha = 0.53
beta = -2.72
gamma = 0.74
hg = 0.024
hr = 0.019

Phi = 1.618
r_4 = 10
r_2h = r_4/(Phi**3)
r_2s = r_4/Phi
R_rot = r_4/(Phi**2)
b_4 = 0.5
b_6 = b_4/Phi
c = r_4/Phi
d = r_4-r_2h
e = r_4/(Phi**2)-b_6
f = r_4 - R_rot
L_imp = r_2h+c+b_6+e
a = c-b_4
b = r_4-r_2s
delta_x = b_6/3
N_bld = 9
e_bld = 0.1
R_rot = 5
L_ind = 40
beta_4 = -45
beta_2 = -56
beta_2s = -60
r_1 = 20
r_5 = 15
e_bld = 0.25
e_tip = 0.01
e_back = 0.01

'''
SGTB Construction
.parameters:
    input Element(dictionary)

.parameters.manual:
    input Length(list), DO3(list), position(integer), alpha(float),
    beta(float), gamma(float), hg(float), hr(float), Ri(float),
    Rg(float), R0(float), L(float)
    should not be used together with .parameters

.grooves:
    input number of grooves(integer)
    by defauly it is 28

.CAD:
    generates the CAD of SGTB for grooves towards right
    for color input 'color' in the end
    for section view input 'section view' in the end
    for dramatized groove depth for applications such as 3D printing input 'dramatize' in the end

.mirror:
    mirrors the right SGTB to get the left SGTB

.combined:
    should be used together with .mirror
    for saving as stl input 'stl' in the end

.right:
    for saving as stl input 'stl' in the end

.left:
    should be used together with .mirror
    for saving as stl input 'stl' in the end
'''

DesignTurbocompressor = Helper()

Element = DesignTurbocompressor.importpickle('Element_23_08_19')

DesignSGTB = SGTB()

DesignSGTB.parameters(Element)
# DesignSGTB.parameters_manual(Length,DO3,pos,alpha,beta,gamma,hg,hr)
DesignSGTB.grooves(28)
DesignSGTB.CAD('color')
DesignSGTB.mirror()
SGTBs = DesignSGTB.combined('stl')
SGTB_right = DesignSGTB.right('stl')
SGTB_left = DesignSGTB.left('stl')

show_object(SGTBs, name='SGTBs')
# show_object(SGTB_right, name='SGTB Right')
# show_object(SGTB_left, name='SGTB Left')

'''
Rotor Construction
.parameters:
    input Element(dictionary)

.parameters.manual:
    input Length(list), DI1(list), DI2(list), DI3(list),
    DO1(list), DO2(list), DO3(list)
    should not be used together with .parameters
    element types can be given in the end by specifying as
    elem_type1 = (list), elem_type2 = (list), elem_type3=(list)

.CAD:
    generates the CAD of rotor
    for color input 'color' in the end
    for section view input 'section view' in the end
    for saving as stl input 'stl' in the end
'''

DesignRotor = Rotor()

DesignRotor.parameters(Element)
# DesignRotor.parameters_manual(Length,DI1,DI2,DI3,DO1,DO2,DO3,elem_type1=types1,elem_type2=types2,elem_type3=types3)
Rotor = DesignRotor.CAD('color','stl')

show_object(Rotor, name='Rotor')


'''
Impeller Construction
.parameters_impeller:
    - extracts the parameters from a pickle file
    - input Element(dictionary)

.manualparams_impeller:
    - user defines parameters manually
    - input Element(dictionary), r_4(float), r_2s(float) ,beta_4(float) ,b_4(float) ,r_1(float), r_2h(float),
      r_5(float), e_bld(float), e_tip(float), e_back(float), L_ind(float), beta_2(float), beta_2s(float), N_bld(integer), R_rot(float)
    - for a custom rotor radius input 'manual_rotor' in the end
    - for using the rotor radius from the pickle file input 'auto_rotor' in the end
    - should not be used together with .parameters_impeller

.hub:
    - models the hub
    - for section view input 'section view' in the end

.blades_excel:
    - retrieves the coordinates of the blades from an excel file
    - input excel 'filename'(string)
    - should not be used together with .blades_coords

.blades_coords:
    - extracts the coordinates of the blade curves from a pickle file
    - input Element(dictionary)
    - should not be used together with .blades_excel

.model_blades:
    - models the blades
    - used after .blades_excel
    - input result of .blades_excel

.rotate_blade(self,blade,bladename):
    - patterns the blades
    - input result of .model_blades, 'bladename'(string)

.assemble(self,files,*settings)
    - combines the impeller components in a common assembly and exports
    - input result of .hub, results of .rotate_blade
    - for exporting  as an atl file input 'stl' or 'STL' in the end
    - exports by default as step
'''

Imp = Impeller()
Imp.parameters_impeller(Element)

#to be used to only model the impeller independently
# Imp.manualparams_impeller(Element,r_4,r_2s,beta_4,b_4,r_1,r_2h,r_5,e_bld,e_tip,e_back,L_ind,beta_2,beta_2s,N_bld,R_rot,'auto_rotor')

Hub = Imp.hub()
# show_object(Hub)

# Coords_mainblades = Imp.blades_excel('coordinates_blade_python.xlsx')
# Coords_splitterblades = Imp.blades_excel('coordinates_splitter_python.xlsx')

Coords_mainblades, Coords_splitterblades = Imp.blades_coords(Element)

Mainblade = Imp.model_blades(Coords_mainblades)
Mainblades = Imp.rotate_blade(Mainblade,'Main Blade')

Splitterblade = Imp.model_blades(Coords_splitterblades)
Splitterblades = Imp.rotate_blade(Splitterblade,'Splitter Blade')

Compressor = Imp.assemble((Hub,Mainblades,Splitterblades))

show_object(Compressor, name = 'Compressor')

DesignTurbocompressor.assemble((Rotor,SGTBs,Hub,Mainblades,Splitterblades),'Turbocompressor')

print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')