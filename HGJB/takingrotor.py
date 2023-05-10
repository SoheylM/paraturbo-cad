import cadquery as cq
import os

cwf = os.getcwd().replace("\\", "/")

# Use this in CQ-Editor
Rotor = cq.importers.importStep(cwf + '/Rotor.stp')

# Use this in VS Code
# Rotor = cq.importers.importStep(cwf + '/HGJB/Rotor.stp')

Rotor=Rotor.transformed(rotate=(0,0,0))

show_object(Rotor)