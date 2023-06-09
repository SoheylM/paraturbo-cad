import pickle
import numpy as np
import cadquery as cq
import os
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions



cwf = os.getcwd().replace("\\", "/")

# Use this in CQ-Editor
Rotor = cq.importers.importStep(cwf + '/RotorExportForComp.step')
CATIA = cq.importers.importStep(cwf + '/SC1-1010-REV1-ROTOR.step')
# Use this in VS Code
# Rotor = cq.importers.importStep(cwf + '/HGJB/Rotor.stp')


Rotor = Rotor.rotate((0,0,0),(0,1,0),90)




#Rotor = Rotor.rotate((0,0,0),(1,0,0),-270)

#cq.exporters.export(Rotor, "RotorExportForComp.step")

show_object(Rotor)

show_object(CATIA)

