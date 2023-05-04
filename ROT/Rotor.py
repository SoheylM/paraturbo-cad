import cadquery as cq
import numpy as np
import os

class Rotor():
    def __init__(self):
        self.cwf = os.getcwd()
        self.method = 'good'

    def parameters(self,Element):
        if type(Element) == dict:
            # Check for information exist in dictionary
            if all(x in Element.keys() for x in ['Laenge','DI1','DI2','DI3','DA1','DA2','DA3','elem_type1','elem_type2','elem_type3']):
                    if len(Element['Laenge']) == len(Element['DI1']) == len(Element['DI2'])  == len(Element['DI3']) == len(Element['DA1']) == len(Element['DA2']) == len(Element['DA3'])\
                        == len(Element['elem_type1']) == len(Element['elem_type2']) == len(Element['elem_type3']):
                        # Taking variables from Element
                        self.length = 1000*np.array(Element['Laenge'])
                        self.DI1 = 1000*np.array(Element['DI1'])
                        self.DI2 = 1000*np.array(Element['DI2'])
                        self.DI3 = 1000*np.array(Element['DI3'])
                        self.DO1 = 1000*np.array(Element['DA1'])
                        self.DO2 = 1000*np.array(Element['DA2'])
                        self.DO3 = 1000*np.array(Element['DA3'])
                        self.elem_type1 = Element['elem_type1']
                        self.elem_type2 = Element['elem_type2']
                        self.elem_type3 = Element['elem_type3']
                        self.elem_type = [self.elem_type1, self.elem_type2, self.elem_type3]
                    else:
                        raise ValueError('Rotor.parameters: Size of the needed dictionary values are not equal.')
            else:
                raise KeyError('Rotor.parameters: Element dictionary does not include all the needed keys.')
        else:
            raise TypeError('Rotor.parameters: Element type is not dictionary.')
        
    def parameters_manual(self,Length,DI1,DI2,DI3,DO1,DO2,DO3,**elemtypes):
        elemtypes_check = True
        if elemtypes:
            if 'elem_type1' in elemtypes and 'elem_type2' in elemtypes and 'elem_type3' in elemtypes:
                if len(elemtypes['elem_type1']) == len(elemtypes['elem_type2']) == len(elemtypes['elem_type3']) == len(Length) == len(DI1) == len(DI2)  == len(DI3) == len(DO1) == len(DO2) == len(DO3):
                    if type(elemtypes['elem_type1']) == list and type(elemtypes['elem_type2']) == list and  type(elemtypes['elem_type3']) == list:
                        for types in elemtypes['elem_type1']:
                            if types in ['COMP1','MAG','ROT','PLUG']:
                                pass
                            else:
                                elemtypes_check = False
                                raise ValueError('Rotor.parameters_manual: Wrong element type is given.')
                        for types in elemtypes['elem_type2']:
                            if types in ['COMP1','MAG','ROT','PLUG']:
                                pass
                            else:
                                elemtypes_check = False
                                raise ValueError('Rotor.parameters_manual: Wrong element type is given.')
                        for types in elemtypes['elem_type3']:
                            if types in ['COMP1','MAG','ROT','PLUG']:
                                pass
                            else:
                                elemtypes_check = False
                                raise ValueError('Rotor.parameters_manual: Wrong element type is given.')
                        if elemtypes_check == True:
                            self.elem_type1 = elemtypes['elem_type1']
                            self.elem_type2 = elemtypes['elem_type2']
                            self.elem_type3 = elemtypes['elem_type3']
                            self.elem_type = [self.elem_type1, self.elem_type2, self.elem_type3]
                    else:
                        raise TypeError('Rotor.parameters_manual: The type of the element type variables are not suitable.')
                else:
                    raise ValueError('Rotor.parameters_manual: Size of the given element types are not equal to each other or not consistent with the length of the other parameters.')
            else:
                raise KeyError('Rotor.parameters_manual: Element types are not given correctly.')
        else:
            self.method = 'manual'

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
                    raise ValueError('Rotor.parameters_manual: Size of the given lists are not equal.')
        else:
            raise TypeError('Rotor.parameters_manual: The type of the given variables are not suitable.')
        
    def CAD(self,*settings):
        # Checking for section view
        if 'section view' in settings:
            sectionview = True
        else:
            sectionview = False

        # Creating dictionaries for layers
        layer1 = {}
        layer2 = {}
        layer3 = {}

        # Creating a count value for further use in union
        count_PLUG = True
        count_ROT = True
        count_MAG = True

        wp = cq.Workplane('XY')

        for i in range(len(self.length)):
            # Constructing 3rd layer
            if self.DO3[i] != self.DI3[i]:
                cylinder3 = wp.cylinder(self.length[i],self.DO3[i]/2,
                direct=(0,0,1),angle=360,centered=(True,True,False),
                combine=False,clean=True).faces('>Z').hole(self.DI3[i],depth=self.length[i],clean=True)
                if sectionview == True:
                    cylinder3 = cylinder3.rect(80,40,(-40,0)).cutThruAll()
                layer3[i] = cylinder3
            else:
                layer3[i] = 0
            # Constructing 2nd layer
            if self.DO2[i] != self.DI2[i]:
                cylinder2 = wp.cylinder(self.length[i],self.DO2[i]/2,
                direct=(0,0,1),angle=360,centered=(True,True,False),
                combine=False,clean=True).faces('>Z').hole(self.DI2[i],depth=self.length[i],clean=True)
                if sectionview == True:
                    cylinder2 = cylinder2.rect(80,40,(-40,0)).cutThruAll()
                layer2[i] = cylinder2
            else:
                layer2[i] = 0
            # Constructing 1st layer
            if self.DO1[i] != self.DI1[i]:
                cylinder1 = wp.cylinder(self.length[i],self.DO1[i]/2,
                direct=(0,0,1),angle=360,centered=(True,True,False),
                combine=False,clean=True).faces('>Z').hole(self.DI1[i],depth=self.length[i],clean=True)
                if sectionview == True:
                    cylinder1 = cylinder1.rect(80,40,(-40,0)).cutThruAll()
                layer1[i] = cylinder1
            else:
                layer1[i] = 0 
            
            # Shifting workplane
            if self.DO1[i] != self.DI1[i]:
                wp = layer1[i].faces('>Z').workplane().center(0,0)
            elif self.DO2[i] != self.DI2[i]:
                wp = layer2[i].faces('>Z').workplane().center(0,0)
            else:
                wp = layer3[i].faces('>Z').workplane().center(0,0)

        layers = [layer1, layer2, layer3]

        if self.method == 'manual':
            assembly = cq.Assembly(name = 'Rotor')
            color = ('red3','green4','blue3','gray50')
            for i in range(0,len(self.length)):
                for k in range(0,3):
                    if layers[k][i] != 0:
                        if 'color' in settings:
                            assembly.add(layers[k][i], name = 'layer'+str(i+1)+'_'+str(k+1), color=cq.Color(color[k]))
                        else:
                            assembly.add(layers[k][i], name = 'layer'+str(i+1)+'_'+str(k+1), color=cq.Color(color[3]))

            assembly.save(self.cwf  + '/STEP/Rotor.step')
        
        else:
            # Uniting cylinders according to their type
            for i in range(0,len(self.length)):
                for k in range(0,3):
                    if layers[k][i] != 0:
                        if self.elem_type[k][i] == 'MAG':
                            if count_MAG == True:
                                MAG = layers[k][i]
                                count_MAG = False
                            else:
                                MAG = MAG.union(layers[k][i])
                        if self.elem_type[k][i] == 'ROT':
                            if count_ROT == True:
                                ROT = layers[k][i]
                                count_ROT = False
                            else:
                                ROT = ROT.union(layers[k][i])
                        if self.elem_type[k][i] == 'PLUG':
                            if count_PLUG == True:
                                PLUG = layers[k][i]
                                count_PLUG = False
                            else:
                                PLUG = PLUG.union(layers[k][i])

            if 'color' in settings:
                assembly = cq.Assembly(name = 'Turbocompressor Rotor')
                assembly.add(ROT, name = 'Rotor', color=cq.Color('green4'))
                assembly.add(PLUG, name = 'Plug', color=cq.Color('blue3'))
                assembly.add(MAG, name = 'Magnet', color=cq.Color('red3'))
            else:
                assembly = cq.Assembly(name = 'Turbocompressor Rotor')
                assembly.add(ROT, name = 'Rotor', color=cq.Color('gray50'))
                assembly.add(PLUG, name = 'Plug', color=cq.Color('gray50'))
                assembly.add(MAG, name = 'Magnet', color=cq.Color('gray50'))

            assembly.save(self.cwf  + '/STEP/Rotor.step')
        
        return assembly