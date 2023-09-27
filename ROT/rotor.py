import cadquery as cq
import numpy as np
import os
import cq_warehouse.extensions

class ROTOR():
    def __init__(self, method_hgjb='Tim', helper_instance=None):
        self.cwf         = os.getcwd().replace("\\", "/")
        self.method      = 'Joseph'
        self.method_hgjb = method_hgjb
        self.helper      = helper_instance

    def parameters(self,Element):
        # Checks for information exists in the dictionary
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
                        self.pos_hgjb2 = Element['sys_pos']['pos_hgjb2']
                        self.alpha_HG = [np.float64(Element['parameters']['hgjb1']['alpha']) , np.float64(Element['parameters']['hgjb2']['alpha'])]
                        self.beta_HG1 = np.float64(Element['parameters']['hgjb1']['beta'])
                        self.beta_HG2 = np.float64(Element['parameters']['hgjb2']['beta'])
                        self.beta_HG = [(self.beta_HG1*np.pi)/180 , (self.beta_HG2*np.pi)/180]
                        self.gamma_HG = [np.float64(Element['parameters']['hgjb1']['gamma']) , np.float64(Element['parameters']['hgjb2']['gamma'])]
                        self.h_gr = [np.float64(Element['parameters']['hgjb1']['hg']) , np.float64(Element['parameters']['hgjb2']['hg'])]
                        self.h_rr = [np.float64(Element['parameters']['hgjb1']['hr']) , np.float64(Element['parameters']['hgjb2']['hr'])]

                        if self.method_hgjb == 'Joseph':
                            # Extract grooves coordinates
                            n_coord = len(Element['parameters']['hgjb2']['x_first_curve'])
                            inner_1 = Element['parameters']['hgjb2']['x_first_curve'][0:n_coord//4]
                            inner_2 = Element['parameters']['hgjb2']['y_first_curve'][0:n_coord//4]
                            inner_3 = Element['parameters']['hgjb2']['z_first_curve'][0:n_coord//4]
                            inner1 = list(zip(inner_1,inner_2,inner_3))
                            self.inner1 = inner1

                            inner_11 = Element['parameters']['hgjb2']['x_first_curve'][3*n_coord//4:]
                            inner_21 = Element['parameters']['hgjb2']['y_first_curve'][3*n_coord//4:]
                            inner_31 = Element['parameters']['hgjb2']['z_first_curve'][3*n_coord//4:]
                            inner11 = list(zip(inner_11,inner_21,inner_31))
                            inner11.reverse()
                            self.inner11 = inner11

                            # Extracting bottom groove surface points
                            P2 = Element['parameters']['hgjb2']['P2_second_surface']
                            self.bottom_groove_surf = P2 #[[cq.Vector(P2[i, j, 0], P2[i, j, 1], P2[i, j, 2]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]

                            # Calculating the shift in the position of the true impeller model
                            self.shift = 0
                            for i in range(len(self.length)):
                                if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                                    self.shift += self.length[i]
                            # Collecting impeller geometry for shift computation 
                            # Rotor radius
                            self.R_rot = (Element['parameters']['comp1']['Rrot'])*1000
                            # Components positions
                            self.pos_comp1 = Element['sys_pos']['pos_comp1']
                            # Tip radius (7mm-35mm)
                            self.r_4 = round(Element['parameters']['comp1']['r4'],12)*1000
                            # Tip width (0.1mm-10.5mm)
                            self.b_4 = (Element['parameters']['comp1']['b4'])*1000
                            # Inlet hub radius (0.7mm-10.5mm)
                            self.r_2h = (Element['parameters']['comp1']['r2h'])*1000
                            # Inlet shoulder radius (0.84mm-24.5mm)
                            self.r_2s = (Element['parameters']['comp1']['r2s'])*1000
                            # Inducer inlet radius (0.84mm-35mm)
                            self.r_1 = (Element['parameters']['comp1']['r1'])*1000
                            # Diffuser exit radius (7mm-52.5mm)
                            self.r_5 = (Element['parameters']['comp1']['r5'])*1000
                            # Blade thickness (0.1mm-0.5mm)
                            self.e_bld = (Element['parameters']['comp1']['Blade_e'])*1000
                            # Inducer length (7mm-140mm)
                            self.L_ind = (Element['parameters']['comp1']['L_ind'])*1000
                            # Tip clearance (0.001mm-0.158mm)
                            self.e_tip = (Element['parameters']['comp1']['Clearance'])*1000
                            # Backface clearance (0.007mm-5.25mm)
                            self.e_back = (Element['parameters']['comp1']['Backface'])*1000
                            # Exit blade angle (-45deg-0deg)
                            self.beta_4 = Element['parameters']['comp1']['beta4']
                            # Defining calculated geometrical parameters
                            self.phi = 1.618
                            self.b_6 = self.b_4/self.phi
                            self.a = (self.r_4/self.phi)-self.b_4
                            self.b = self.r_4-self.r_2s
                            self.c = self.r_4/self.phi
                            self.d = self.r_4-self.r_2h
                            self.e =(self.r_4/(self.phi**2))-self.b_6
                            self.f = self.r_4-self.R_rot
                            self.L_imp = self.r_2h+self.c+self.b_6+self.e
                            # Offset due to impeller : comp_offset = Limp-shift
                            self.comp_offset = self.L_imp - self.shift #1.6355 # need to grab the actual values if necessary

                            
                    else:
                        raise ValueError('ROTOR.parameters: Size of the needed dictionary values are not equal.')
            else:
                raise KeyError('ROTOR.parameters: Element dictionary does not include all the needed keys.')
        else:
            raise TypeError('ROTOR.parameters: Element type is not a dictionary.')
        
    def parameters_manual(self,Length,DI1,DI2,DI3,DO1,DO2,DO3,pos_HGJB1,pos_HGJB2,alpha_HGJB,beta_HGJB,gamma_HGJB,hg_HGJB,hr_HGJB,**elemtypes):
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
                                raise ValueError('ROTOR.parameters_manual: Wrong element type is given.')
                        for types in elemtypes['elem_type2']:
                            if types in ['COMP1','MAG','ROT','PLUG']:
                                pass
                            else:
                                elemtypes_check = False
                                raise ValueError('ROTOR.parameters_manual: Wrong element type is given.')
                        for types in elemtypes['elem_type3']:
                            if types in ['COMP1','MAG','ROT','PLUG']:
                                pass
                            else:
                                elemtypes_check = False
                                raise ValueError('ROTOR.parameters_manual: Wrong element type is given.')
                        if elemtypes_check == True:
                            self.elem_type1 = elemtypes['elem_type1']
                            self.elem_type2 = elemtypes['elem_type2']
                            self.elem_type3 = elemtypes['elem_type3']
                            self.elem_type = [self.elem_type1, self.elem_type2, self.elem_type3]
                    else:
                        raise TypeError('ROTOR.parameters_manual: The type of the element type variables are not suitable.')
                else:
                    raise ValueError('ROTOR.parameters_manual: Size of the given element types are not equal to each other or not consistent with the length of the other parameters.')
            else:
                raise KeyError('ROTOR.parameters_manual: Element types are not given correctly.')
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
                    self.pos_hgjb1 = np.int64(pos_HGJB1)
                    self.pos_hgjb2 = np.int64(pos_HGJB2)
                    self.alpha_HG = [np.float64(alpha_HGJB)]
                    self.beta_HG = np.float64(beta_HGJB)
                    self.beta_HG = [(np.pi*self.beta_HG)/180]
                    self.gamma_HG = [np.float64(gamma_HGJB)]
                    self.h_gr = [np.float64(hg_HGJB)]
                    self.h_rr = [np.float64(hr_HGJB)]
                else:
                    raise ValueError('ROTOR.parameters_manual: Size of the given lists are not equal.')
        else:
            raise TypeError('ROTOR.parameters_manual: The type of the given variables are not suitable.')
        
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

        self.layers = [layer1, layer2, layer3]

        # Enters if the element types are not given correctly or not given
        if self.method == 'manual':
            for i in range(0,len(self.length)):
                for k in range(0,3):
                    if self.layers[k][i] != 0:
                        if count_ROT == True:
                            ROT = self.layers[k][i]
                            count_ROT = False
                        else:
                            ROT = ROT.union(self.layers[k][i])
            
            self.ROT = ROT

        # Enters if the elements types are given correctly or a dictionary is used
        else:
            # Uniting cylinders according to their type
            for i in range(0,len(self.length)):
                for k in range(0,3):
                    if self.layers[k][i] != 0:
                        if self.elem_type[k][i] == 'MAG':
                            if count_MAG == True:
                                MAG = self.layers[k][i]
                                count_MAG = False
                            else:
                                MAG = MAG.union(self.layers[k][i])
                        if self.elem_type[k][i] == 'ROT':
                            if count_ROT == True:
                                ROT = self.layers[k][i]
                                count_ROT = False
                            else:
                                ROT = ROT.union(self.layers[k][i])
                        if self.elem_type[k][i] == 'PLUG':
                            if count_PLUG == True:
                                PLUG = self.layers[k][i]
                                count_PLUG = False
                            else:
                                PLUG = PLUG.union(self.layers[k][i])

            self.PLUG = PLUG
            self.MAG = MAG
            self.ROT = ROT

        return ROT
    
    def HGJB(self):
        self.N_HG = 28                 # Number of grooves generally between 26 - 30
        self.D = self.DO3[self.pos_hgjb1] #16                    # On drawing [mm]
        # Length of HGJB on drawing [mm]
        self.L = [self.length[self.pos_hgjb1], self.length[self.pos_hgjb2]]
        # Value for CAD
        self.L_land=[self.L[0]-(self.gamma_HG[0]*self.L[0]) , self.L[1]-(self.gamma_HG[1]*self.L[1])]
        # Oversized length for safety generally between 0.6 - 1
        self.L= [self.L[0] + 0.8 , self.L[1] + 0.8]                                                                          
        self.Spiral_step = [np.pi*self.D*np.tan(self.beta_HG[0]), np.pi*self.D*np.tan(self.beta_HG[1])]
        self.Spiral_height = [self.L[0]/2 , self.L[1]/2]
        self.a_HG = [(np.pi*self.D*self.alpha_HG[0])/self.N_HG , (np.pi*self.D*self.alpha_HG[1])/self.N_HG]
        self.a_HG_plus_b_HG = [self.a_HG[0]/self.alpha_HG[0] , self.a_HG[1]/self.alpha_HG[1]]
        # Diametral clearance given in micrometers
        self.h_rr_tot = [self.h_rr[0]*2 , self.h_rr[1]*2]

        self.dist1 = 0
        for d in range(self.pos_hgjb1):
            self.dist1=self.dist1+self.length[d]
            
        # Find distance to first HGJB
        self.dist2 = 0
        for d in range(self.pos_hgjb2):
            self.dist2=self.dist2+self.length[d]

        # Find distance to center of first HGJB
        self.DistCenter1=self.dist1+self.L[0]/2

        # Find distance to center of second HGJB
        self.DistCenter2=self.dist2+self.L[1]/2

        # Define separation angle between grooves
        self.sepang=360/self.N_HG

        # Adding discretization
        self.n_parall = 10
        self.eps_perc = 0.005

        # Length between parallelogram verticals
        self.LenBetwVert = [(self.L[0]/2 - self.L_land[0]/2)/self.n_parall , (self.L[1]/2 - self.L_land[1]/2)/self.n_parall]

        # Angle between parallelogram diagonal and reference horizontal
        self.Betaprime = [abs(self.beta_HG[0]) - np.pi/2 , abs(self.beta_HG[1]) - np.pi/2]

        # Gap between horizontal leaving from parallelogram lowest 
        # Corner and parallelogram vertical

        self.gap = [self.LenBetwVert[0]*np.tan(self.Betaprime[0]) , self.LenBetwVert[1]*np.tan(self.Betaprime[1])]
        self.gap_spiral = [(self.gap[0]*self.Spiral_step[0])/self.LenBetwVert[0] , (self.gap[1]*self.Spiral_step[1])/self.LenBetwVert[1]]

        if self.method_hgjb == 'Joseph':
            # Distance to translate hgjb1
            self.dist_hack0 = -(self.DistCenter1+self.comp_offset-self.L_land[0]/2) #-50
            # Distance to increase L_land : set to zero
            self.dist_hack1 = 0 #1
            # Distance to translate hgjb2
            self.dist_hack2 = self.DistCenter2 -(self.DistCenter1) #60
            # Step angle in degrees
            self.angle_step = self.sepang #360/self.N_HG 

    def HGJB_CAD(self,rotor):

        if self.method_hgjb == 'Tim':
            rotor = rotor.rotate((0,0,0),(1,0,0),270)

            # Create cylinder parameters
            self.CylLen1=self.length[self.pos_hgjb1]
            self.CylRadOut1 = 1.01*self.DO3[self.pos_hgjb1]/2
            
            # Projection cylinder
            self.cylinder1 = cq.Solid.makeCylinder(
                self.CylRadOut1, 2*self.CylLen1+2, pnt=cq.Vector(0,1*self.L[0],0), dir=cq.Vector(0,-1,0)
                )
            
            # Direction of projection
            self.projection_direction = cq.Vector(0,0,-1)
            
            # Calculate turn angle in radians
            self.radang = self.gap[0]/self.CylRadOut1
            # Convert to degrees
            self.rotang = self.radang*180/np.pi
            
            self.parallelograms1 = []
            self.parallelograms2 = []
            self.parallelograms_projected1 = []
            self.parallelograms_projected2 = []
            
            for k in range(0,2):
                for i in range(self.n_parall):
                    self.parallelogram = (
                        cq.Workplane("YX", origin = ((self.gap[k]+self.a_HG[k])/2, -self.L[k]/2+self.LenBetwVert[k] +i*self.LenBetwVert[k],0))
                        # Drawing parallelogram in local coordinates
                        .lineTo(-(self.LenBetwVert[k])*(1+self.eps_perc),-self.gap[k])
                        .lineTo(-(self.LenBetwVert[k])*(1+self.eps_perc),-self.gap[k]-self.a_HG[k])
                        .lineTo(0,-self.a_HG[k])
                        .close()
                        .extrude(1)
                        .faces("<Z")
                        .val()
                        )
                    if k == 0: 
                        self.parallelograms1.append(self.parallelogram)
                    if k == 1:
                        self.parallelograms2.append(self.parallelogram)
                
            # Project onto cylinder
            self.para_seg = []
            for k in range(0,2):
                for i in range(self.n_parall):
                    if i == 0:
                        self.cylinder1 = self.cylinder1.rotate((0,0,0), (0,1,0),self.rotang)

                    if k == 0:
                        self.parallelogram_projected = self.parallelograms1[i].projectToShape(self.cylinder1, self.projection_direction)
                        self.parallelograms_projected1.append(self.parallelogram_projected)
                    if k == 1:
                        self.parallelogram_projected = self.parallelograms2[i].projectToShape(self.cylinder1, self.projection_direction)
                        self.parallelograms_projected2.append(self.parallelogram_projected)
                    
                    if i == self.n_parall-1:
                        pass
                    else:
                        self.cylinder1 = self.cylinder1.rotate((0,0,0), (0,1,0), self.rotang)
            
                if k == 0:
                    self.para_solid1 = cq.Compound.makeCompound(
                        [f1.thicken(self.h_gr[k]*1000*9, cq.Vector(0,0,1)) for f1 in self.parallelograms_projected1[0]]
                        )
                if k == 1:
                    self.para_solid2 = cq.Compound.makeCompound(
                        [f2.thicken(self.h_gr[k]*1000*9, cq.Vector(0,0,1)) for f2 in self.parallelograms_projected1[0]]  # should be 2 but causes problem
                        )
            
                if k == 0:
                    self.para_seg.append(self.para_solid1)
                    for j in range(self.n_parall-1):
                        self.para_seg_tr1 = self.para_seg[0].transformed((0, -(j+1)*self.rotang, 0), (0, (j+1)*self.LenBetwVert[0],0))
                        self.para_solid1 = self.para_solid1.fuse(self.para_seg_tr1)
                if k == 1:
                    self.para_seg.append(self.para_solid2)
                    for j in range(self.n_parall-1):
                        self.para_seg_tr2 = self.para_seg[1].transformed((0, -(j+1)*self.rotang, 0), (0, (j+1)*self.LenBetwVert[1],0))
                        self.para_solid2 = self.para_solid2.fuse(self.para_seg_tr2)
                    
            # Mirror about the origin plane
            self.para_solid_m = [self.para_solid1.mirror('XZ'), self.para_solid2.mirror('XZ')]
            
            # Creating actual HGJB1 and HGJB2 groove objects, near is closer to impeller than far
            self.groove1near = {}
            self.groove1near[0] = self.para_solid1.transformed((0,0, 0), (0, self.DistCenter1, 0))

            self.groove1far = {}
            self.groove1far[0] = self.para_solid_m[0].transformed((0,0, 0), (0, self.DistCenter1, 0))

            self.groove2near = {}
            self.groove2near[0] = self.para_solid2.transformed((0,0, 0), (0, self.DistCenter2, 0))

            self.groove2far = {}
            self.groove2far[0] = self.para_solid_m[1].transformed((0,0, 0), (0, self.DistCenter2, 0))

            # Actually cutting the rotor
            for i in range(0,2):
                self.groove1near[i+1] = self.groove1near[i].transformed ((0 ,self.sepang ,0))
                rotor = rotor.cut(self.groove1near[i+1])
        
                self.groove1far[i+1] = self.groove1far[i].transformed ((0 ,self.sepang ,0))
                rotor = rotor.cut(self.groove1far[i+1])
        
                self.groove2near[i+1] = self.groove2near[i].transformed ((0 ,self.sepang ,0))
                rotor = rotor.cut(self.groove2near[i+1])
        
                self.groove2far[i+1] = self.groove2far[i].transformed ((0 ,self.sepang ,0))
                rotor = rotor.cut(self.groove2far[i+1])

            rotor = rotor.rotate((0,0,0),(1,0,0),-270)

            self.ROT = rotor
        
        elif self.method_hgjb == 'Joseph':
            print('method_hgjb ==', self.method_hgjb)
            # Load the required libraries
            from math import pi, tan
            from cqmore import Workplane
            from cqmore.matrix import translation, rotationZ
            from cqmore.matrix import mirror, translation

            # Initialize dictionaries
            profile = {}
            profile1 = {}

            usage={}
            usagemm={}
            surface={}

            usa={}
            mir={}
            mir1={}
            surf={}

            usa2={}
            mir2={}
            mir12={}
            surf2={}

            surfacemm={}

            mirrored_pts={}
            mirrored_pts1={}

            trans={}
            trans1={}
            tra={}
            tra1={}

            # Extract previously computed variables
            inner1     = self.inner1
            inner11    = self.inner11
            dist_hack0 = self.dist_hack0
            dist_hack1 = self.dist_hack1
            dist_hack2 = self.dist_hack2
            angle_step = self.angle_step

            # Translate the rotor
            rotor = rotor.translate((0,0,dist_hack0))

            ## NEW METHOD ##
            # Define half bearing with pooints on the surface:
            P2 = self.bottom_groove_surf

            # Convert P2 to a flat list of tuples
            #hgjb1_half1 = P2.reshape(-1, 3).tolist()
            hgjb1_half1 = [tuple(point) for point in P2.reshape(-1, 3).tolist()]


            # Define transformations
            t2 = translation((0, 0, dist_hack2)) # translate to hgjb2
            mirror_operation = mirror((0, 0, 1)) # mirror half bearing

            # Apply transformations
            hgjb1_half2 = mirror_operation.transformAll(hgjb1_half1) # mirror groove at hgjb1 (2 half bearings)
            hgjb2_half1 = t2.transformAll(hgjb1_half1) 
            hgjb2_half2 = t2.transformAll(hgjb1_half2)

            # Initialize dictionaries
            dic_hgjb1_half1 = {}
            dic_hgjb1_half2 = {}
            dic_hgjb2_half1 = {}
            dic_hgjb2_half2 = {}

            surf_hgjb1_half1 = {}
            surf_hgjb1_half2 = {}
            surf_hgjb2_half1 = {}
            surf_hgjb2_half2 = {}

            # Assign first values
            dic_hgjb1_half1[0] = hgjb1_half1
            dic_hgjb1_half2[0] = hgjb1_half2
            dic_hgjb2_half1[0] = hgjb2_half1
            dic_hgjb2_half2[0] = hgjb2_half2

            # Define rotation matrix
            m = rotationZ(angle_step)

            # Now you can use these with splineApproxSurface
            thickness = 0.1  # example value for thickness, adjust as needed

            # Convert transformed points to 2D list of Vector objects            
            points_hgjb1_half1 = [[cq.Vector(*hgjb1_half1[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]
            points_hgjb1_half2 = [[cq.Vector(*hgjb1_half2[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]
            points_hgjb2_half1 = [[cq.Vector(*hgjb2_half1[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]
            points_hgjb2_half2 = [[cq.Vector(*hgjb2_half2[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]

            # Perform loop to cover N_HG=28 grooves
            for i in range(0,self.N_HG):

                dic_hgjb1_half1[i+1] = m.transformAll(dic_hgjb1_half1[i])
                hgjb1_half1_tmp      = dic_hgjb1_half1[i]
                points_hgjb1_half1 = [[cq.Vector(*hgjb1_half1_tmp[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]
                surf_hgjb1_half1[i+1]  = Workplane().splineApproxSurface(points_hgjb1_half1,-thickness,clean=True,combine=False)
                #print('dic_hgjb1_half1 i =', i)

                dic_hgjb1_half2[i+1] = m.transformAll(dic_hgjb1_half2[i])
                hgjb1_half2_tmp      = dic_hgjb1_half2[i]
                points_hgjb1_half2 = [[cq.Vector(*hgjb1_half2_tmp[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]
                surf_hgjb1_half2[i+1]  = Workplane().splineApproxSurface(points_hgjb1_half2,thickness,clean=True,combine=False)
                #print('dic_hgjb1_half2 i =', i)

                dic_hgjb2_half1[i+1] = m.transformAll(dic_hgjb2_half1[i])
                hgjb2_half1_tmp      = dic_hgjb2_half1[i]
                points_hgjb2_half1 = [[cq.Vector(*hgjb2_half1_tmp[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]
                surf_hgjb2_half1[i+1]  = Workplane().splineApproxSurface(points_hgjb2_half1,-thickness,clean=True,combine=False)
                #print('dic_hgjb2_half1 i =', i)

                dic_hgjb2_half2[i+1] = m.transformAll(dic_hgjb2_half2[i])
                hgjb2_half2_tmp      = dic_hgjb2_half2[i]
                points_hgjb2_half2 = [[cq.Vector(*hgjb2_half2_tmp[i * P2.shape[1] + j]) for j in range(P2.shape[1])] for i in range(P2.shape[0])]
                surf_hgjb2_half2[i+1]  = Workplane().splineApproxSurface(points_hgjb2_half2,thickness,clean=True,combine=False)
                #print('dic_hgjb2_half2 i =', i)

            texts = Workplane()
            for i in range(0,self.N_HG): 
                texts.add(surf_hgjb1_half1[i+1])
                texts.add(surf_hgjb1_half2[i+1])
                texts.add(surf_hgjb2_half1[i+1])
                texts.add(surf_hgjb2_half2[i+1])
                #print('texts add i=', i)


            # Subract the grooves from the rotor
            rotor = rotor - texts

            # Translate back the rotor
            rotor = rotor.translate((0,0,-dist_hack0))

            # Testing export
            cq.exporters.export(rotor,self.cwf  + '/STEP/Rotor_Joseph.step', opt={"write_pcurves": False, "precision_mode": -1})

        #self.ROT = rotor
        
    def assemble(self,*settings):
       # Enters if the element types are not given correctly or not given
        if self.method == 'manual':
            ROT = cq.importers.importStep(self.cwf + '/STEP/Rotor_Joseph.step') #self.ROT
            assembly = cq.Assembly(name = 'Turbocompressor Rotor')
            assembly.add(ROT, name = 'Rotor', color=cq.Color('gray50')) 

        else:
            PLUG = self.PLUG
            ROT = cq.importers.importStep(self.cwf + '/STEP/Rotor_Joseph.step')
            MAG = self.MAG
            # Adds parts to assembly with our without color
            if 'color' in settings:
                print('color wanted on rotor.')
                assembly = cq.Assembly(name = 'Turbocompressor Rotor')
                assembly.add(ROT, name   = 'Rotor', color=cq.Color('green4'))
                assembly.add(PLUG, name  = 'Plug', color=cq.Color('blue3'))
                assembly.add(MAG, name   = 'Magnet', color=cq.Color('red3'))
            else:
                assembly = cq.Assembly(name = 'Turbocompressor Rotor')
                assembly.add(ROT, name   = 'Rotor', color=cq.Color('gray50'))
                print('assembly ROT')
                assembly.add(PLUG, name  = 'Plug', color=cq.Color('gray50'))
                print('assembly PLUG')
                assembly.add(MAG, name   = 'Magnet', color=cq.Color('gray50'))
                print('assembly MAG')

            # Saves as stl if given in arguments
            #if 'stl' or 'STL' in settings:
            #    cq.exporters.export(ROT, self.cwf + '/STL/Rotor.stl')
            #    cq.exporters.export(PLUG, self.cwf + '/STL/Plug.stl')
            #    cq.exporters.export(MAG, self.cwf + '/STL/Magnet.stl')
        
        # Saves as step
        assembly.save(self.cwf  + '/STEP/Rotor.step')
        #cq.exporters.export(assembly, self.cwf  + '/STEP/Rotor.step', opt={"write_pcurves": False, "precision_mode": -1})
        print('assembly saved')

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Rotor', self.cwf + '/STL/Rotor')

        return assembly