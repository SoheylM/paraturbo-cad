import cadquery as cq
import numpy as np
import os
import cq_warehouse.extensions

class Rotor():
    def __init__(self):
        self.cwf = os.getcwd().replace("\\", "/")
        self.method = 'Joseph'

    def parameters(self,Element):
        # Checks for information exist in dictionary
        # Checks for are all needed keys are in Element, are all values have same lenght, is Element a dictionary
        # If not raises an error
        if type(Element) == dict:
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
                        # Taking variables for HGJB
                        self.pos_hgjb1 = Element['sys_pos']['pos_hgjb1']
                        self.pos_hgjb2 =Element['sys_pos']['pos_hgjb2']
                        self.alpha_HG = [0.68]
                        self.beta_HG = -135
                        self.beta_HG = [(self.beta_HG*np.pi)/180]
                        self.gamma_HG = [0.89]
                        self.h_gr = [16/1000]
                        self.h_rr = [9/1000]
                        # self.alpha_HG = [np.float64(Element['parameters']['hgjb1']['alpha']) , np.float64(Element['parameters']['hgjb2']['alpha'])]
                        # self.beta_HG1 = np.float64(Element['parameters']['hgjb1']['beta'])
                        # self.beta_HG2 = np.float64(Element['parameters']['hgjb2']['beta'])
                        # self.beta_HG = [(self.beta_HG1*np.pi)/180 , (self.beta_HG2*np.pi)/180]
                        # self.gamma_HG = [np.float64(Element['parameters']['hgjb1']['gamma']) , np.float64(Element['parameters']['hgjb2']['gamma'])]
                        # self.h_gr = [1000*np.float64(Element['parameters']['hgjb1']['hg']) , 1000*np.float64(Element['parameters']['hgjb2']['hg'])]
                        # self.h_rr = [1000*np.float64(Element['parameters']['hgjb1']['hr']) , 1000*np.float64(Element['parameters']['hgjb2']['hr'])]
                    else:
                        raise ValueError('Rotor.parameters: Size of the needed dictionary values are not equal.')
            else:
                raise KeyError('Rotor.parameters: Element dictionary does not include all the needed keys.')
        else:
            raise TypeError('Rotor.parameters: Element type is not dictionary.')
        
    def parameters_manual(self,Length,DI1,DI2,DI3,DO1,DO2,DO3,**elemtypes):
        # Checks for element types if provided
        # Checks for the length and type of the given element types and the inputs for the element types
        # If not raises an error
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

        # Checks for information exist in given input
        # Checks for are all values have same lenght and type
        # If not raises an error
        if type(Length) == list and type(DO3) == list and type(DO2) == list and type(DO1) == list and type(DI3) == list and type(DI2) == list and type(DI1) == list: 
                if len(Length) == len(DI1) == len(DI2)  == len(DI3) == len(DO1) == len(DO2) == len(DO3):
                    # Taking variables from user
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
        # Checking for section view in arguments
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

        # Creating the workplane
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

        # Enters if the element types are not given correctly or not given
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

            # Later convert this to raise Error and combine with HGJB code!!!           
            print('Rotor.CAD: Rotor could not be provided due to missing element types.')

        # Enters if the elements types are given correctly or a dictionary is used
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

            # Adds parts to assembly with our without color
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

        # Saves as stl if given in arguments
        if 'stl' or 'STL' in settings:
            cq.exporters.export(ROT, self.cwf + '/STL/Rotor.stl')
            cq.exporters.export(PLUG, self.cwf + '/STL/Plug.stl')
            cq.exporters.export(MAG, self.cwf + '/STL/Magnet.stl')
            # cq.exporters.export(assembly.toCompound(), self.cwf + '/STL/Turbocompressor Rotor.stl')
        
        # Saves as step
        # cq.exporters.export(ROT, self.cwf + '/HGJB/Rotor.stp', cq.exporters.ExportTypes.STEP)
        assembly.save(self.cwf  + '/STEP/Rotor.step')

        return assembly, ROT
    
    def HGJB(self,rotor):
        self.N_HG = 28 #number of grooves generally between 26 - 30
        self.D = 16 #on drawing [mm]
        self.L = self.length[self.pos_hgjb1]#28 #length of HGJB on drawing [mm]
        self.L_land=self.L-(self.gamma_HG[0]*self.L) #Value for CAD
        self.L=self.L+0.8 #oversized length for safety generally between 0.6 - 1
        self.Spiral_step = np.pi*self.D*np.tan(self.beta_HG[0])
        self.Spiral_height = self.L/2
        self.a_HG = (np.pi*self.D*self.alpha_HG[0])/self.N_HG #mm
        self.a_HG_plus_b_HG = self.a_HG/self.alpha_HG[0] #mm
        self.h_rr_tot = self.h_rr[0]*2 #diametral clearance given in micrometers
        #end code adapted from Christophe's Matlab

        self.dist1 = 0
        for d in range(self.pos_hgjb1):
            self.dist1=self.dist1+self.length[d]
            
        #find distance to first HGJB
        self.dist2 = 0
        for d in range(self.pos_hgjb2):
            self.dist2=self.dist2+self.length[d]

        #find distance to center of first HGJB
        self.DistCenter1=self.dist1+self.L/2

        #find distance to center of second HGJB
        self.DistCenter2=self.dist2+self.L/2

        #define separation angle between grooves
        self.sepang=360/self.N_HG

        #length between parallelogram verticals
        self.LenBetwVert = self.L/2 - self.L_land/2

        #Angle between parallelogram diagonal and reference horizontal
        self.Betaprime = abs(self.beta_HG[0]) - 90

        #gap between horizontal leaving from parallelogram lowest 
        #corner and parallelogram vertical

        self.gap = self.LenBetwVert*np.tan(self.Betaprime)

        #create removal cylinder1
        self.CylLen1 = self.length[self.pos_hgjb1]
        self.CylRadOut1= self.DO3[self.pos_hgjb1]/2

        #create removal cylinder2
        self.CylLen2 = self.length[self.pos_hgjb2]
        self.CylRadOut2= self.DO3[self.pos_hgjb2]/2
        self.rot = rotor

    def HGJB_CAD(self):
        self.rot = self.rot.rotate((0,0,0),(1,0,0),270)
        #find distance to first HGJB
        
        #first cylinder to be projected onto
        cylinder1 = cq.Solid.makeCylinder(
             self.CylRadOut1, self.CylLen1+2, pnt=cq.Vector(0, self.DistCenter1+self.L/2+1, 0), dir=cq.Vector(0, -1, 0)
        )
        #first removal cylinder
        removalcylinder1 = cq.Solid.makeCylinder(
              self.CylRadOut1*2, self.CylLen1, pnt=cq.Vector(0, self.DistCenter1+self.L/2, -self.DO3[self.pos_hgjb1]*0.8), dir=cq.Vector(0, -1, 0)
        )

        #second cylinder to be projected onto
        cylinder2 = cq.Solid.makeCylinder(
             self.CylRadOut2, self.CylLen2+2, pnt=cq.Vector(0, self.DistCenter2+self.L/2+1, 0), dir=cq.Vector(0, -1, 0)
        )
        #second removal cylinder
        removalcylinder2 = cq.Solid.makeCylinder(
              self.CylRadOut2*2, self.CylLen2, pnt=cq.Vector(0, self.DistCenter2+self.L/2, -self.DO3[self.pos_hgjb1]*0.8), dir=cq.Vector(0, -1, 0)
        )

        #direction of projection 
        projection_direction = cq.Vector(0, 0, 1)

        #global coordinates of origin of parallelogram1
        # yp1 = DistCenter1+L/2
        xp1 = self.LenBetwVert*np.tan(self.Betaprime)/2 
        yp2 = -self.LenBetwVert*np.tan(self.Betaprime)
        xp2 = -self.LenBetwVert


        parallelogram1a = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter1+self.L/2, -2*self.CylRadOut1))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(-self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(-self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        parallelogram2a = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter1-self.L/2, -2*self.CylRadOut1))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        parallelogram1b = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter2+self.L/2, -2*self.CylRadOut2))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(-self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(-self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        parallelogram2b = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter2-self.L/2, -2*self.CylRadOut2))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        #fist HGJB
        #project first parallelogram onto cylinder
        parallelogram1a_projected = parallelogram1a.projectToShape(cylinder1, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram1a_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr[0]) for f in parallelogram1a_projected]
          )
        parallelogram1a_solids = parallelogram1a_solids.cut(removalcylinder1)

        #project second parallelogram onto cylinder
        parallelogram2a_projected = parallelogram2a.projectToShape(cylinder1, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram2a_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr[0]) for f in parallelogram2a_projected]
          )
        parallelogram2a_solids = parallelogram2a_solids.cut(removalcylinder1)


        #second HGJB
        #project first parallelogram onto cylinder
        parallelogram1b_projected = parallelogram1b.projectToShape(cylinder2, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram1b_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr[0]) for f in parallelogram1b_projected]
          )
        parallelogram1b_solids = parallelogram1b_solids.cut(removalcylinder2)

        #project second parallelogram onto cylinder
        parallelogram2b_projected = parallelogram2b.projectToShape(cylinder2, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram2b_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr[0]) for f in parallelogram2b_projected]
          )
        parallelogram2b_solids = parallelogram2b_solids.cut(removalcylinder2)

        for i in range(self.N_HG):
            self.rot = self.rot.cut(parallelogram1a_solids)
            self.rot = self.rot.cut(parallelogram2a_solids)
            self.rot = self.rot.cut(parallelogram1b_solids)
            self.rot = self.rot.cut(parallelogram2b_solids)
            self.rot = self.rot.rotate((0,0,0),(0,1,0),self.sepang)

        self.rot = self.rot.rotate((0,0,0),(1,0,0),-270)
        
        return self.rot