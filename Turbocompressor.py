import cadquery as cq
import numpy as np
import time
import os
cwf = os.getcwd()

import sys
sys.path.append(cwf  + '/SGTB')
sys.path.append(cwf  + '/ROT')
sys.path.append(cwf  + '/COMP')

from SGTB import *
from Rotor import *
from Impeller import *
from Helper import *

'''
IMPORTANT
SAVE CODE BEFORE RUNNING OTHERWISE .STEP FILE WILL BE WRONGLY GENERATED
'''

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
pos = 6
alpha = 0.53
beta = -2.72
gamma = 0.74
hg = 0.024
hr = 0.019

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
.combined:
    should be used together with .mirror
.left:
    should be used together with .mirror
'''

DesignTurbocompressor = Helper()

Element = DesignTurbocompressor.importpickle('Element_23_04_22')

DesignSGTB = SGTB()

DesignSGTB.parameters(Element)
# DesignSGTB.parameters_manual(Length,DO3,pos,alpha,beta,gamma,hg,hr)
DesignSGTB.grooves(28)
DesignSGTB.CAD('color')
DesignSGTB.mirror()
SGTBs = DesignSGTB.combined()
# SGTB_right = DesignSGTB.right()
# SGTB_left = DesignSGTB.left()

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
.CAD:
    generates the CAD of rotor
    for color input 'color' in the end
    for section view input 'section view' in the end
'''

DesignRotor = Rotor()

DesignRotor.parameters(Element)
# DesignRotor.parameters_manual(Length,DI1,DI2,DI3,DO1,DO2,DO3)
Rotor = DesignRotor.CAD('color')

show_object(Rotor, name='Rotor')

Imp = Impeller()

Imp.parameters_impeller(Element)

Imp.settings_hub(True,True,False)
Hub = Imp.hub()

#modeling the main blades
Coords_mainblades = Imp.blades_excel('POINT_BLADES1.xls')
Mainblade = Imp.model_blades(Coords_mainblades)
Mainblades = Imp.rotate_blade(Mainblade)

#modeling the splitter blades
Coords_splitterblades = Imp.blades_excel('POINT_BLADES2.xls')
Splitterblade = Imp.model_blades(Coords_splitterblades)
Splitterblades = Imp.rotate_blade(Splitterblade)

show_object(Hub)
show_object(Mainblades)
show_object(Splitterblades)

DesignTurbocompressor.assemble((SGTBs,Rotor,Hub,Mainblades,Splitterblades),'Turbocompressor')

print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')