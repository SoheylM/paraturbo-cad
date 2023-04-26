import cadquery as cq
import numpy as np
import pandas as pd
from cadquery.selectors import AreaNthSelector
import pickle
from cadquery import *
from cq_warehouse import *
from cq_warehouse.extensions import *
import os

class Impeller():
    def __init__(self):
        
        #default impeller variables (in mm) if function is called without parameters from the pickle file

        #tip radius (7mm-35mm)
        self.r_4 = 19

        #inlet shoulder radius (0.84mm-24.5mm)
        self.r_2s = 2

        #tip width (0.1mm-10.5mm)
        self.b_4 = 2

        #inducer inlet radius (0.84mm-35mm)
        self.r_1 = 20

        #inducer hub radius (0.7mm-10.5mm)
        self.r_2h = 3

        # diffuser exit radius (7mm-52.5mm)
        self.r_5 = 15

        #blade thickness (0.1mm-0.5mm)
        self.e_bld = 0.25

        #tip clearance (0.001mm-0.158mm)
        self.e_tip = 0.01

        #backface clearance (0.007mm-5.25mm)
        self.e_back = 0.01

        #inducer length (7mm-140mm)
        self.L_ind = 40

        #exit blade angle (-45deg-0deg)
        self.beta_4 = -45

        #preventing the twistExtrude() function from failing
        if self.beta_4 ==0:
            self.beta_4 = 0.01

        #inlet blade angle (-45deg)
        self.beta_2 = -56

        #inlet blade angle shroud (-56deg)
        self.beta_2s = -60

        #number of blades (5blades-11blades)
        self.N_bld = 9

        # radius of rotor
        self.R_rot = 5

        #defining geometrical parameters
        self.phi = 1.618
        self.b_6 = self.b_4/self.phi
        self.a = (self.r_4/self.phi)-self.b_4
        self.b = self.r_4-self.r_2s
        self.c = self.r_4/self.phi
        self.d = self.r_4-self.r_2h
        self.e =(self.r_4/(self.phi**2))-self.b_6
        self.f = self.r_4-self.R_rot
        self.L_imp = self.r_2h+self.c+self.b_6+self.e

        self.top_aux=True
        self.bottom_aux=True
        self.top_impeller=True
        self.bottom_impeller=True
        self.top_hub = True
        self.bottom_hub = True

    def importpickle(self,filename):
        self.cwf = os.path.dirname(os.path.abspath(__file__))
        file = open(self.cwf  + '/' + filename, 'rb')
        Element = pickle.load(file)
        file.close

        return Element
    
    def parameters(self,Element):
        
        #################################################
        #extracting element parameters from pickle in mm#
        #################################################

        #the coordinate system origin is the center of reference

        #diamaters and heights of cylindrical representation
        self.Laenge = [i * 1000 for i in Element['Laenge']]
        self.DI1 = [i * 1000 for i in Element['DI1']]
        self.DI2 = [i * 1000 for i in Element['DI2']]
        self.DI3 = [i * 1000 for i in Element['DI3']]
        self.DA1 = [i * 1000 for i in Element['DA1']]
        self.DA2 = [i * 1000 for i in Element['DA2']]
        self.DA3 = [i * 1000 for i in Element['DA3']]

        #tip radius (7mm-35mm)
        self.r_4 = round(Element['parameters']['comp1']['r4'],12)*1000

        #tip width (0.1mm-10.5mm)
        self.b_4 = (Element['parameters']['comp1']['b4'][0])*1000

        #inlet hub radius (0.7mm-10.5mm)
        self.r_2h = (Element['parameters']['comp1']['r2h'])*1000

        #inlet shoulder radius (0.84mm-24.5mm)
        self.r_2s = (Element['parameters']['comp1']['r2s'])*1000

        #inducer inlet radius (0.84mm-35mm)
        self.r_1 = (Element['parameters']['comp1']['r1'])*1000

        # diffuser exit radius (7mm-52.5mm)
        self.r_5 = (Element['parameters']['comp1']['r5'])*1000

        #blade thickness (0.1mm-0.5mm)
        self.e_bld = (Element['parameters']['comp1']['Blade_e'])*1000

        #inducer length (7mm-140mm)
        self.L_ind = (Element['parameters']['comp1']['L_ind'])*1000

        #tip clearance (0.001mm-0.158mm)
        self.e_tip = (Element['parameters']['comp1']['Clearance'])*1000

        #backface clearance (0.007mm-5.25mm)
        self.e_back = (Element['parameters']['comp1']['Backface'])*1000

        #exit blade angle (-45deg-0deg)
        self.beta_4 = Element['parameters']['comp1']['beta4']

        if self.beta_4 ==0:
            self.beta_4 = 0.01

        #number of blades (5blades-11blades)
        self.N_bld = int(Element['parameters']['comp1']['N_blades'])

        # radius of rotor
        self.R_rot = (Element['parameters']['comp1']['R_ROT'])*1000

        #inlet blade angle (-45deg) (constant)
        beta_2 = -56

        #inlet blade angle shroud (-56deg) (constant)
        beta_2s = -60
            
        #defining geometrical parameters
        self.phi = 1.618
        self.b_6 = self.b_4/self.phi
        self.a = (self.r_4/self.phi)-self.b_4
        self.b = self.r_4-self.r_2s
        self.c = self.r_4/self.phi
        self.d = self.r_4-self.r_2h
        self.e =(self.r_4/(self.phi**2))-self.b_6
        self.f = self.r_4-self.R_rot
        self.L_imp = self.r_2h+self.c+self.b_6+self.e

    def settings(self,top_a,bot_a,top_i,bot_i,top_h,bot_h):
        self.top_aux=top_a
        self.bottom_aux=bot_a
        self.top_impeller=top_i
        self.bottom_impeller=bot_i
        self.top_hub = top_h
        self.bottom_hub = bot_h

    def hub(self):
        #######################################
        #modeling the hub from the pickle file#
        #######################################

        #controls the smoothness of the curve
        # step_hub=int(L_imp/b_6)*10
        self.step_hub=100

        #initializing the x and y coordinate arrays
        self.x_hub = np.linspace(0,self.L_imp,self.step_hub)
        self.y_hub=[]
        self.y_hub.insert(0,0)

        #dividing the hub into 4 parts: hub_1,hub_2,hub_3,hub_4
        #looping over the 4 parts to update the x and y cooridinates
        for i in range(0,len(self.x_hub)):

        #hub_1
            if (self.x_hub[i]>=0)&(self.x_hub[i]<=self.r_2h):
                self.y_hub[i] = eval('((self.r_2h**2)-((self.x_hub[i]-self.r_2h)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])

        #hub_2
            if (self.x_hub[i]>=self.r_2h)&(self.x_hub[i]<=self.r_2h+self.c):
                self.y_hub[i] = eval('self.r_4-((self.d**2)-((self.d/self.c)**2)*((self.x_hub[i]-self.r_2h)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])

        #hub_3
            if (self.x_hub[i]>=self.r_2h+self.c)&(self.x_hub[i]<=self.r_2h+self.c+self.b_6):
                self.y_hub[i]=self.r_4
                self.y_hub.insert(i,self.y_hub[i])

        #hub_4
            if (self.x_hub[i]>=self.r_2h+self.c+self.b_6)&(self.x_hub[i]<=self.L_imp):

                self.y_hub[i] = eval('self.r_4-((self.f**2)-((self.f/self.e)**2)*((self.x_hub[i]-self.L_imp)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])

        #grouping the x and y coordinates into tuples
        self.coords_hub = list(zip(self.x_hub,self.y_hub))

        # skteching half the hub continuously
        self.sketch_hub = (cq.Workplane("XY")
                    .polyline(self.coords_hub,includeCurrent=False)
                    .vLineTo(0)
                    .close()
                    )

        #revolving the sketch about the x axis
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                .translate((0,0,0))
                .split(keepTop=self.top_hub,keepBottom=self.bottom_hub)
                )
        
        return hub