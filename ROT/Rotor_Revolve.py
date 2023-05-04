#importing required libraries
import cadquery as cq
import numpy as np
import pandas as pd
from cadquery.selectors import AreaNthSelector
from cadquery import *
from cq_warehouse import *
from cq_warehouse.extensions import *
import os

#defining a class to model the impeller
class Rotor_Revolve():

    #defining the class constructor
    def __init__(self):
        self.top_layer1=True
        self.bottom_layer1=True
        self.top_layer2=True
        self.bottom_layer2=True
        self.top_layer3=True
        self.bottom_layer3=True
        self.auto_Rrot = True
        self.cwf = os.getcwd()

    #defining a method to extract rotor parameters from the pickle file in mm
    def parameters_rotor(self,Element):

        #diamaters and heights of cylindrical representation
        self.Laenge = [i * 1000 for i in Element['Laenge']]
        self.DI1 = [i * 1000 for i in Element['DI1']]
        self.DI2 = [i * 1000 for i in Element['DI2']]
        self.DI3 = [i * 1000 for i in Element['DI3']]
        self.DA1 = [i * 1000 for i in Element['DA1']]
        self.DA2 = [i * 1000 for i in Element['DA2']]
        self.DA3 = [i * 1000 for i in Element['DA3']]

        #component types by layer
        self.elem_type1 = Element['elem_type1']
        self.elem_type2 = Element['elem_type2']
        self.elem_type3 = Element['elem_type3']

        #components positions
        self.pos_comp1 = Element['sys_pos']['pos_comp1']

    #defining a method to model the rotor
    def model_rotor(self):

        #initialising disctionaries to store geometrical shapes per function
        layer1={}
        layer2={}
        layer3={}

        #initialising the position at the center of each element along the x axis
        hor_pos=[]
        shift=0
        hor_pos.insert(0,self.Laenge[0]/2)

        #modeling every element using a loop and the revolve function
        for i in range(len(self.Laenge)):

            #not modeling the cylindrical representation of the impeller
            if self.elem_type1[i]!='COMP1' or self.elem_type2[i]!='COMP1' or self.elem_type3[i]!='COMP1':
                #modeling layer 1
                #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                    if self.DI1[i]!=self.DA1[i]:
                        layer1[i] = (cq.Workplane("XY")
                                .moveTo(hor_pos[i],self.DI1[i]/2+(self.DA1[i]-self.DI1[i])/4)
                                .rect(self.Laenge[i],(self.DA1[i]-self.DI1[i])/2)
                                .revolve(360,(0,0,0),(1,0,0))
                                .split(keepTop=self.top_layer1,keepBottom=self.bottom_layer1)
                                .rotate((0,0,0),(0,1,0),-90))
                            
                #modeling layer 2
                #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                    if self.DI2[i]!=self.DA2[i]:
            
                        layer2[i] = (cq.Workplane("XY")
                                    .moveTo(hor_pos[i],self.DI2[i]/2+(self.DA2[i]-self.DI2[i])/4)
                                    .rect(self.Laenge[i],(self.DA2[i]-self.DI2[i])/2)
                                    .revolve(360,(0,0,0),(1,0,0))
                                    .split(keepTop=self.top_layer2,keepBottom=self.bottom_layer2)
                                    .rotate((0,0,0),(0,1,0),-90))
            
                #modeling layer 3
                #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                    if self.DI3[i]!=self.DA3[i]:
            
                        layer3[i] = (cq.Workplane("XY")
                                        .moveTo(hor_pos[i],self.DI3[i]/2+(self.DA3[i]-self.DI2[i])/4)
                                        .rect(self.Laenge[i],(self.DA3[i]-self.DI3[i])/2)
                                        .revolve(360,(0,0,0),(1,0,0))
                                        .split(keepTop=self.top_layer2,keepBottom=self.bottom_layer3)
                                        .rotate((0,0,0),(0,1,0),-90))
                                        

        #updating the position at the center of each element along the x axis
            if i<len(self.Laenge)-1:
                hor_pos.insert(i+1,hor_pos[i]+self.Laenge[i]/2+self.Laenge[i+1]/2)

        #asembling the three layers of the rotor
        assembly = cq.Assembly()
        for j in range(len(self.Laenge)):

            if j in layer1.keys():
                assembly.add(layer1[j],color=cq.Color("red"))

            if j in layer2.keys():
                assembly.add(layer2[j],color=cq.Color("green"))

            if j in layer3.keys():
                assembly.add(layer3[j],color=cq.Color("blue"))

        return assembly
    
    #defining a method to change the rotor modeling settings
    def settings_rotor(self,t_layer1,b_layer1,t_layer2,b_layer2,t_layer3,b_layer3):
        #booleans to create a section view
        self.top_layer1=t_layer1
        self.bottom_layer1=b_layer1
        self.top_layer2=t_layer2
        self.bottom_layer2=b_layer2
        self.top_layer3=t_layer3
        self.bottom_layer3=b_layer3