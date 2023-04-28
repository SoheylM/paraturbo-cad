
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
class Impeller():

    #defining the class constructor
    def __init__(self):
        self.top_hub = True
        self.bottom_hub = True
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

    #defining a method to extract the impeller parameters from the pickle file in mm
    def parameters_impeller(self,Element):     

        #tip radius (7mm-35mm)
        self.r_4 = round(Element['parameters']['comp1']['r4'],12)*1000

        #tip width (0.1mm-10.5mm)
        self.b_4 = (Element['parameters']['comp1']['b4'])*1000

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

        '''
        if self.beta_4 == 0:
            self.beta_4 = 0.01
        '''

        #number of blades (5blades-11blades)
        self.N_bld = int(Element['parameters']['comp1']['N_blades'])

        #inlet blade angle (constant)
        self.beta_2 = -56

        #inlet blade angle shroud (constant)
        self.beta_2s = -60

        #component types by layer
        self.elem_type1 = Element['elem_type1']
        self.elem_type2 = Element['elem_type2']
        self.elem_type3 = Element['elem_type3']

        #rotor radius
        self.R_rot = (Element['parameters']['comp1']['Rrot'])*1000

        #defining calculated geometrical parameters
        self.phi = 1.618
        self.b_6 = self.b_4/self.phi
        self.a = (self.r_4/self.phi)-self.b_4
        self.b = self.r_4-self.r_2s
        self.c = self.r_4/self.phi
        self.d = self.r_4-self.r_2h
        self.e =(self.r_4/(self.phi**2))-self.b_6
        self.f = self.r_4-self.R_rot
        self.L_imp = self.r_2h+self.c+self.b_6+self.e

        #components positions
        self.pos_comp1 = Element['sys_pos']['pos_comp1']

    #defining a method to manually define impeller variables (in mm) by the user
    def manualparams_impeller(self,r_4,r_2s,beta_4,b_4,r_1,r_2h,r_5,e_bld,e_tip,e_back,L_ind,beta_2,beta_2s,N_bld,R_rot):

        self.r_4 = r_4 
        self.r_2s = r_2s
        self.beta_4 = beta_4 
        self.b_4 = b_4 
        self.rr_1 = r_1 
        self.r_2h = r_2h 
        self.r_5 = r_5 
        self.e_bld = e_bld 
        self.ce_tip = e_tip 
        self.e_back = e_back 
        self.L_ind = L_ind 
        self.beta_2 = beta_2 
        self.beta_2s = beta_2s 
        self.N_bld = N_bld 
        self.R_rot = R_rot

        if self.auto_Rrot==True:
            #automatically defining the radius of the rotor from the pickle file
            for k in range(len(self.Laenge)):
                if self.elem_type1[k]!='COMP1' or self.elem_type2[k]!='COMP1' or self.elem_type3[k]!='COMP1':
                            break
                self.R_rot = self.DA3[k]/2

        #defining calculated geometrical parameters
        self.phi = 1.618
        self.b_6 = self.b_4/self.phi
        self.a = (self.r_4/self.phi)-self.b_4
        self.b = self.r_4-self.r_2s
        self.c = self.r_4/self.phi
        self.d = self.r_4-self.r_2h
        self.e =(self.r_4/(self.phi**2))-self.b_6
        self.f = self.r_4-self.R_rot
        self.L_imp = self.r_2h+self.c+self.b_6+self.e

        return

    # defining a method to change the hub modeling settings
    def settings_hub(self,top_h,bot_h):

        #boolean parameters to show a section of the hub
        self.top_hub = top_h
        self.bottom_hub = bot_h

    # defining a method to model the hub
    def hub(self):

        r_2h = self.r_2h
        r_4 = self.r_4
        b_4 = self.b_4
        
        phi = 1.618
        b_6 = b_4/phi
        c = r_4/phi
        e =(r_4/(phi**2))-b_6
        L_imp = r_2h+c+b_6+e

        # tolerance to control the smoothness of the hub profile
        self.step_hub=100

        #initializing the x and y coordinate arrays
        self.x_hub = np.linspace(0,L_imp,self.step_hub)
        self.y_hub=[]
        self.y_hub.insert(0,0)

        #dividing the hub into 4 parts: hub_1,hub_2,hub_3,hub_4
        #looping over the 4 parts to update the x and y coordinates
        for i in range(0,len(self.x_hub)):

        #hub_1
            if (self.x_hub[i]>=0)&(self.x_hub[i]<=r_2h):
                self.y_hub[i] = eval('((self.r_2h**2)-((self.x_hub[i]-self.r_2h)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])

        #hub_2
            if (self.x_hub[i]>=r_2h)&(self.x_hub[i]<=r_2h+c):
                self.y_hub[i] = eval('self.r_4-((self.d**2)-((self.d/self.c)**2)*((self.x_hub[i]-self.r_2h)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])

        #hub_3
            if (self.x_hub[i]>=r_2h+c)&(self.x_hub[i]<=r_2h+c+b_6):
                self.y_hub[i]=r_4
                self.y_hub.insert(i,self.y_hub[i])

        #hub_4
            if (self.x_hub[i]>=r_2h+c+b_6)&(self.x_hub[i]<=L_imp):

                self.y_hub[i] = eval('self.r_4-((self.f**2)-((self.f/self.e)**2)*((self.x_hub[i]-self.L_imp)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])

        #grouping the x and y coordinates into tuples
        self.coords_hub = list(zip(self.x_hub,self.y_hub))

        # skteching half the hub continuously
        self.sketch_hub = (cq.Workplane("XY")
                    .polyline(self.coords_hub,includeCurrent=False)
                    .vLineTo(0)
                    .close())

        # calculating the shift in the position of the true impeller model
        shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                shift += self.Laenge[i]

        #revolving the sketch about the rotation axis
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                .translate((-abs(self.L_imp-shift),0,0))
                .split(keepTop=self.top_hub,keepBottom=self.bottom_hub))
    
        assembly = cq.Assembly()
        assembly.add(hub,color=cq.Color('red'))
        return assembly
    
    # defining a method to retrieve the coordinates of the blades from the excel file
    def blades_excel(self,filename):

        #storing lists of non-planar blade curves
        points_blade = {}
        coords_blade = {}
        
        #reading the excel file with blade points
        excelfile = self.cwf  + '/COMP/' + filename
        print(excelfile)
        df_blade = pd.read_excel(io=excelfile,header=None)
        
        #setting index markers to indicate start and end of storing values
        start = None
        end = None
        
        #looping to iterate over all the rows
        for index, row in df_blade.iterrows():
        
        #updating index markers every time 'StartCurve' or 'EndCurve' is detected
            if 'StartCurve' in row.values:
                start = index+ 1 #skip the header row
            elif 'EndCurve' in row.values:
                end = index
        
        #storing curve points
                if len(points_blade) == 0:
                    points_blade[0] = df_blade.iloc[start:end,:]
                else:
                    points_blade[len(points_blade)] = df_blade.iloc[start:end,:]
                start = None
                end = None
        
        #grouping the curve x,y,z coordinates in tuples
        for i in range(len(points_blade)):
            coords_blade[i] = list(zip(points_blade[i].iloc[:,0].astype(float),\
                                            points_blade[i].iloc[:,1].astype(float),\
                                            points_blade[i].iloc[:,2].astype(float)))
        
        return coords_blade
    
    #defining a method to model the blades
    def model_blades(self,coords):
                
        #Adding the first point of the curve to the end of its point array to close it
        coords[0].append(coords[0][0])
        coords[len(coords)-1].append(coords[len(coords)-1][0])
        
        #initializing lists to store the solids
        blade_solid={}
        
        #sketching the first blade curve as an edge
        blade_edge0 = cq.Edge.makeSpline([cq.Vector(p) for p in coords[0]][0:-1]).close()
        
        #creating a closed surface within the first blade curve edge
        blade_face0 = cq.Face.makeNSidedSurface([blade_edge0],[])
        
        #sketching the last blade curve as an edge
        blade_edge_last = cq.Edge.makeSpline([cq.Vector(p) for p in coords[len(coords)-1]][0:-1]).close()
        
        #creating a closed surface within the last blade curve edge
        blade_face_last = cq.Face.makeNSidedSurface([blade_edge_last],[])
        
        #lofting the closed surfaces to produce the overall blade shape
        blade_lofted = cq.Solid.makeLoft([blade_edge0,blade_edge_last],True)
        
        #shelling the faces and edges
        blade_shell = cq.Shell.makeShell([blade_face0,blade_face_last,blade_lofted]).fix()
        
        #solidifying the produced shell and rotating
        blade_solid[0] = cq.Solid.makeSolid(blade_shell).rotate((0,0,0),(0,1,0),90).translate((-13,0,0))

        return blade_solid
    
    #defining a method to pattern the blades
    def rotate_blade(self,blade):

        assembly = cq.Assembly(blade[0],color=cq.Color('red'))

        #looping to model the blade pattern
        for i in range(0,self.N_bld-1):
        
            #rotating about the x axis by the corresponding angle
            blade[i+1] = (blade[i].transformed((360/self.N_bld,0,0)))
            assembly.add(blade[i+1],color=cq.Color('red'))

        return assembly
    
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
                                .split(keepTop=self.top_layer1,keepBottom=self.bottom_layer1)                                )
                            
                #modeling layer 2
                #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                    if self.DI2[i]!=self.DA2[i]:
            
                        layer2[i] = (cq.Workplane("XY")
                                    .moveTo(hor_pos[i],self.DI2[i]/2+(self.DA2[i]-self.DI2[i])/4)
                                    .rect(self.Laenge[i],(self.DA2[i]-self.DI2[i])/2)
                                    .revolve(360,(0,0,0),(1,0,0))
                                    .split(keepTop=self.top_layer2,keepBottom=self.bottom_layer2)                                    )
            
                #modeling layer 3
                #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                    if self.DI3[i]!=self.DA3[i]:
            
                        layer3[i] = (cq.Workplane("XY")
                                        .moveTo(hor_pos[i],self.DI3[i]/2+(self.DA3[i]-self.DI2[i])/4)
                                        .rect(self.Laenge[i],(self.DA3[i]-self.DI3[i])/2)
                                        .revolve(360,(0,0,0),(1,0,0))
                                        .split(keepTop=self.top_layer2,keepBottom=self.bottom_layer3)                                        )

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
    def settings_rotor(self,t_layer1,b_layer1,t_layer2,b_layer2,t_layer3,b_layer3,bool_Rrot):
        #manually defining the radius of the rotor from the pickle file
        self.top_layer1=t_layer1
        self.bottom_layer1=b_layer1
        self.top_layer2=t_layer2
        self.bottom_layer2=b_layer2
        self.top_layer3=t_layer3
        self.bottom_layer3=b_layer3
        #automatically defining the radius of the rotor from the pickle file
        self.auto_Rrot = bool_Rrot

