import cadquery as cq
import numpy as np
import os
import pickle

class HELPER():
    def __init__(self):
        self.cwf = os.getcwd().replace("\\", "/")

    def importpickle(self,filename):
        # Importing pickle file and closing the file
        file = open(self.cwf  + '/ELEMENT/' + filename + '.pickle', 'rb')
        Element = pickle.load(file)
        file.close

        return Element
    
    def assemble(self,files,filename):
        # Assembling the parts that are given and saving as step
        assembly = cq.Assembly(name=filename)
        for i in range(0,len(files)):
            assembly.add(files[i],name='Subassembly '+str(i+1))

        assembly.save(self.cwf  + '/STEP/' + filename + '.step')