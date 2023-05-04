import cadquery as cq
import numpy as np
import os
import pickle

class Helper():
    def __init__(self):
        self.cwf = os.getcwd()

    def importpickle(self,filename):
        file = open(self.cwf  + '/' + filename + '.pickle', 'rb')
        Element = pickle.load(file)
        file.close

        return Element
    
    def assemble(self,files,filename):
        assembly = cq.Assembly(name=filename)
        for i in range(0,len(files)):
            assembly.add(files[i],name='Subassembly '+str(i+1))

        assembly.save(self.cwf  + '/STEP/' + filename + '.step')