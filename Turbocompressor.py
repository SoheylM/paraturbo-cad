import cadquery as cq
import numpy as np
import time
import os
cwf = os.path.dirname(os.path.abspath(__file__))

import sys
sys.path.append(cwf  + '/SGTB')
sys.path.append(cwf  + '/ROT')

from SGTB import *
from Rotor import *
from GeneralFuncs import *

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
    beta(float), gamma(float), hg(float), hr(float)
    should not be used together with .parameters
.settings:
    input color(boolean) & sectionview(boolean)
    by default it is True,False
.grooves:
    input number of grooves(integer)
    by defauly it is 28
.CAD:
    generates the CAD of SGTB for grooves towards right
.combined:
    should be used together with .mirror & .position
.left:
    should be used together with .mirror
'''

Gen = GeneralFuncs()

Element = Gen.importpickle('Element_23_04_22.pickle')

SGTB = SGTB()

SGTB.parameters(Element)
# SGTB.parameters_manual(Length,DO3,pos,alpha,beta,gamma,hg,hr)
SGTB.settings(True,False)
SGTB.grooves(28)
SGTB.CAD()
SGTB.mirror()
SGTB.position()
my_SGTBs = SGTB.combined()
# my_SGTB_right = SGTB.right()
# my_SGTB_left = SGTB.left()

# show_object(my_SGTBs, name='SGTBs')
# show_object(my_SGTB_right, name='SGTB Right')
# show_object(my_SGTB_left, name='SGTB Left')

'''
Rotor Construction
.parameters:
    input Element(dictionary)
.parameters.manual:
    input Length(list), DI1(list), DI2(list), DI3(list),
    DO1(list), DO2(list), DO3(list)
    should not be used together with .parameters
.settings:
    input color(boolean) & sectionview(boolean)
    by default it is True,False
.CAD:
    generates the CAD of rotor
'''

Rotor = Rotor()

Rotor.parameters(Element)
# Rotor.parameters_manual(Length,DI1,DI2,DI3,DO1,DO2,DO3)
Rotor.settings(True,False)
my_Rotor = Rotor.CAD()

# show_object(my_Rotor, name='Rotor')

Gen.assemble((my_SGTBs,my_Rotor))

print('Time: ' + str(np.round((time.time()-t0),2)) + ' seconds')