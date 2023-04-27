import cadquery as cq
import numpy as np
import os
import pickle

class Helper():
    def __init__(self):
        self.cwf = os.getcwd()

    def importpickle(self,filename):
        file = open(self.cwf  + '/' + filename, 'rb')
        Element = pickle.load(file)
        file.close

        return Element
    
    def assemble(self,files):
        assembly = cq.Assembly()
        for i in range(0,len(files)):
            assembly.add(files[i])

        assembly.save(self.cwf  + '/Turbocompressor.step')