import cadquery as cq
import numpy as np
import os
import pickle
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

class Rotor():
    def __init__(self):
        self.cwf = os.getcwd().replace("\\", "/")
        # Use this in CQ-Editor
        
    def imp(self): 
        return cq.importers.importStep(self.cwf + '/Rotor.stp')
        
        
# cwf = os.getcwd().replace("\\", "/")
# # Use this in CQ-Editor
# Rotor = cq.importers.importStep(cwf + '/Rotor.stp')

t = Rotor().imp()

show_object(t)