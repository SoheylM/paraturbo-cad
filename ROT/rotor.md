# Rotor Code

import cadquery as cq
import numpy as np
import os
import pickle

class Rotor():
    def __init__(self):
        self.color = True
        self.sectionview = False
        self.cwf = os.path.dirname(os.path.abspath(__file__))

    def importpickle(self,filename):
        file = open(self.cwf  + '/' + filename, 'rb')
        Element = pickle.load(file)
        file.close

        return Element
    
    def parameters(self,Element):
        if type(Element) == dict:
            # Check for information exist in dictionary
            if all(x in Element.keys() for x in ['Laenge','DI1','DI2','DI3','DA1','DA2','DA3']):
                    if len(Element['Laenge']) == len(Element['DI1']) == len(Element['DI2'])  == len(Element['DI3']) == len(Element['DA1']) == len(Element['DA2']) == len(Element['DA3']):
                        # Taking variables from Element
                        self.length = 1000*np.array(Element['Laenge'])
                        self.DI1 = 1000*np.array(Element['DI1'])
                        self.DI2 = 1000*np.array(Element['DI2'])
                        self.DI3 = 1000*np.array(Element['DI3'])
                        self.DO1 = 1000*np.array(Element['DA1'])
                        self.DO2 = 1000*np.array(Element['DA2'])
                        self.DO3 = 1000*np.array(Element['DA3'])
                    else:
                        print('Rotor.parameters: Size of the needed dictionary values are not equal.')
                        return
            else:
                print('Rotor.parameters: Element dictionary does not include all the needed keys.')
                return
        else:
            print('Rotor.parameters: Element type is not dictionary.')
            return 
        
    def parameters_manual(self,Length,DI1,DI2,DI3,DO1,DO2,DO3):
        if type(Length) == list and type(DO3) == list and type(DO2) == list and type(DO1) == list and type(DI3) == list and type(DI2) == list and type(DI1) == list: 
                if len(Length) == len(DI1) == len(DI2)  == len(DI3) == len(DO1) == len(DO2) == len(DO3):
                    self.length = np.array(Length)
                    self.DI1 = np.array(DI1)
                    self.DI2 = np.array(DI2)
                    self.DI3 = np.array(DI3)
                    self.DO1 = np.array(DO1)
                    self.DO2 = np.array(DO2)
                    self.DO3 = np.array(DO3)
                else:
                    print('Rotor.parameters_manual: Size of the given list values are not equal.')
                    return
        else:
            print('Rotor.parameters_manual: The type of the given variables are not suitable.')
            return
        
    def settings(self,color,sectionview):
        self.color = color
        self.sectionview = sectionview

    def CAD(self):
        # Creating an empty cylinder dictionary
        cylinders = {}

        # Defining workplane
        start = cq.Workplane('XY')

        # Creating a fuction to make a three layered cylinder
        def threelayercylinder(di1,di2,di3,do1,do2,do3,l,wp,sw):
            if di3 != do3:
                layer3 = wp.cylinder(l,do3/2,
                direct=(0,0,1),angle=360,centered=(True,True,False),
                combine=False,clean=True).faces('>Z').hole(di3,depth=l,clean=True)
                if sw == True:
                    layer3 = layer3.rect(80,40,(-40,0)).cutThruAll()
            else:
                layer3 = wp

            if di2 != do2:
                layer2 = wp.cylinder(l,do2/2,
                direct=(0,0,1),angle=360,centered=(True,True,False),
                combine=False,clean=True).faces('>Z').hole(di2,depth=l,clean=True)
                if sw == True:
                    layer2 = layer2.rect(80,40,(-40,0)).cutThruAll()
            else:
                layer2 = wp

            if di1 != do1:
                layer1 = wp.cylinder(l,do1/2,
                direct=(0,0,1),angle=360,centered=(True,True,False),
                combine=False,clean=True).faces('>Z').hole(di1,depth=l,clean=True)
                if sw == True:
                    layer1 = layer1.rect(80,40,(-40,0)).cutThruAll()
            else:
                layer1 = wp

            return (layer3,layer2,layer1)

        # Using the function in loop to create turbocompressor
        color = ('blue3','green4','red3','gray50')
        assembly = cq.Assembly()
        for i in range(0,len(self.length)):
            cylinders['cylinder'+str(i)] = threelayercylinder(self.DI1[i],self.DI2[i],self.DI3[i],self.DO1[i],self.DO2[i],self.DO3[i],self.length[i],start,self.sectionview)
            locals().update(cylinders)
            if self.DO1[i] != self.DI1[i]:
                start = cylinders['cylinder'+str(i)][2].faces('>Z').workplane().center(0,0)
            elif self.DO2[i] != self.DI2[i]:
                start = cylinders['cylinder'+str(i)][1].faces('>Z').workplane().center(0,0)
            else:
                start = cylinders['cylinder'+str(i)][0].faces('>Z').workplane().center(0,0)

            if self.color == False:
                for k in range(0,3):
                    assembly.add(
                        cylinders['cylinder'+str(i)][k],
                        name='cylinder'+str(i)+str(k),color=cq.Color(color[3]))
            elif self.color == True:
                for k in range(0,3):
                    assembly.add(
                        cylinders['cylinder'+str(i)][k],
                        name='cylinder'+str(i)+str(k),color=cq.Color(color[k]))

        assembly.save(self.cwf  + '/Rotor.step')

        return assembly