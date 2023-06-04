import cadquery as cq
import numpy as np
import os
import pickle
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions




class Rotor():
    def __init__(self):
        #open file
        self.file = open('Element_23_08_19.pickle', 'rb')

        # dump info to that file
        self.Element = pickle.load(self.file)

        #close file
        self.file.close()

        self.Laenge = self.Element['Laenge']
        #change list into numpy array
        self.Laenge = np.array(self.Laenge)
        #change units from m to mm to avoid hollow visual effect

        self.Laenge = 1000*self.Laenge 
        #multiplying Laenge streeeetches the part out, so don't
        self.DI1 = 1000*np.array(self.Element['DI1'])
        self.DI2 = 1000*np.array(self.Element['DI2'])
        self.DI3 = 1000*np.array(self.Element['DI3'])
        self.DA1 = 1000*np.array(self.Element['DA1'])
        self.DA2 = 1000*np.array(self.Element['DA2'])
        self.DA3 = 1000*np.array(self.Element['DA3'])
        #print('type DI1', type(DI1))

        self.sys_pos = self.Element['sys_pos']
        self.pos_hgjb1 = self.sys_pos['pos_hgjb1']
        self.pos_hgjb2 =self.sys_pos['pos_hgjb2']
        self.cwf = os.getcwd().replace("\\", "/")
        #Begin code adapted from Christophe's Matlab
        #Parameters 
        self.N_HG = 28 #number of grooves generally between 26 - 30
        self.alpha_HG = self.Element['parameters']['hgjb1']['alpha']#0.68 #given
        self.beta_HG = self.Element['parameters']['hgjb1']['beta'] #-135 #given 
        self.beta_HG = self.beta_HG*pi/180
        self.gamma_HG = self.Element['parameters']['hgjb1']['gamma'] #0.89 #given
        self.h_gr = self.Element['parameters']['hgjb1']['hg'] #16 #groove depth given in micrometers
        self.h_rr = 9 #clearance on radiu given in micrometers
        self.D = 16 #on drawing [mm]
        self.L = self.Laenge[self.pos_hgjb1]#28 #length of HGJB on drawing [mm]
        self.L_land=self.L-(self.gamma_HG*self.L) #Value for CAD
        self.L=self.L+0.8 #oversized length for safety generally between 0.6 - 1
        self.Spiral_step = pi*self.D*tan(self.beta_HG)
        self.Spiral_height = self.L/2
        self.a_HG = (pi*self.D*self.alpha_HG)/self.N_HG #mm
        self.a_HG_plus_b_HG = self.a_HG/self.alpha_HG #mm
        self.h_rr_tot = self.h_rr*2 #diametral clearance given in micrometers
        #end code adapted from Christophe's Matlab
        
        self.dist1 = 0
        for d in range(self.pos_hgjb1):
            self.dist1=self.dist1+self.Laenge[d]
            
        #find distance to first HGJB
        self.dist2 = 0
        for d in range(self.pos_hgjb2):
            self.dist2=self.dist2+self.Laenge[d]
            
        #find distance to center of first HGJB
        self.DistCenter1=self.dist1+self.L/2

        #find distance to center of second HGJB
        self.DistCenter2=self.dist2+self.L/2

        #define separation angle between grooves
        self.sepang=360/self.N_HG

        #length between parallelogram verticals
        self.LenBetwVert = self.L/2 - self.L_land/2

        #Angle between parallelogram diagonal and reference horizontal
        self.Betaprime = abs(self.beta_HG) - pi/2

        #gap between horizontal leaving from parallelogram lowest 
        #corner and parallelogram vertical

        self.gap = self.LenBetwVert*tan(self.Betaprime)
        self.gap_spiral = (self.gap*self.Spiral_step)/self.LenBetwVert

        #create removal cylinder1
        self.CylLen1 = self.Laenge[self.pos_hgjb1]
        self.CylRadOut1= self.DA3[self.pos_hgjb1]/2

        #create removal cylinder2
        self.CylLen2 = self.Laenge[self.pos_hgjb2]
        self.CylRadOut2= self.DA3[self.pos_hgjb2]/2
        
        #adding discretization
        self.n_parall = 10
        self.eps_perc = 0.005
        
        
        self.rot = cq.importers.importStep(self.cwf + '/Rotor.stp')
        
    def create_HGJB(self): 
        
        self.rot = self.rot.rotate((0,0,0),(1,0,0),270)
        #find distance to first HGJB
        

        

        self.rot = self.rot.rotate((0,0,0),(1,0,0),-270)
        
        
        
        return self.rot
        
        


t = Rotor().create_HGJB()


show_object(t)
