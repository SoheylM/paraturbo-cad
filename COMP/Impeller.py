# Importing required libraries
import cadquery as cq
import numpy as np
import pandas as pd
#from cadquery import *
#from cq_warehouse import *
#from cq_warehouse.extensions import *
import os
#from OCP.TopoDS import TopoDS_Shape 
import time

# Defining a class to model the impeller
class IMPELLER():

    # Defining the class constructor
    def __init__(self, helper_instance=None):
        self.auto_Rrot = True
        self.cwf       = os.getcwd().replace("\\", "/")
        self.helper    = helper_instance

    # Defining a method to extract the impeller parameters from the pickle file in mm
    def parameters_impeller(self,Element):     

        # Diamaters and heights of cylindrical representation
        self.Laenge = [i * 1000 for i in Element['Laenge']]
        self.DI1 = [i * 1000 for i in Element['DI1']]
        self.DI2 = [i * 1000 for i in Element['DI2']]
        self.DI3 = [i * 1000 for i in Element['DI3']]
        self.DA1 = [i * 1000 for i in Element['DA1']]
        self.DA2 = [i * 1000 for i in Element['DA2']]
        self.DA3 = [i * 1000 for i in Element['DA3']]

        # Component types by layer
        self.elem_type1 = Element['elem_type1']
        self.elem_type2 = Element['elem_type2']
        self.elem_type3 = Element['elem_type3']

        # Components positions
        self.pos_comp1 = Element['sys_pos']['pos_comp1']

        # Tip radius (7mm-35mm)
        self.r_4 = round(Element['parameters']['comp1']['r4'],12)*1000 # necesary to avoid bugs
        #self.r_4 = (Element['parameters']['comp1']['r4'])*1000

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

        '''
        if self.beta_4 == 0:
            self.beta_4 = 0.01
        '''

        # Number of blades (5blades-11blades)
        self.N_bld = int(Element['parameters']['comp1']['N_blades'])

        # Inlet blade angle (constant)
        self.beta_2 = -56

        # Inlet blade angle shroud (constant)
        self.beta_2s = -60

        # Component types by layer
        self.elem_type1 = Element['elem_type1']
        self.elem_type2 = Element['elem_type2']
        self.elem_type3 = Element['elem_type3']

        # Rotor radius
        self.R_rot = (Element['parameters']['comp1']['Rrot'])*1000

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

        return

    # Defining a method to manually define impeller variables (in mm) by the user
    def manualparams_impeller(self,Element,r_4,r_2s,beta_4,b_4,r_1,r_2h,r_5,e_bld,e_tip,e_back,L_ind,beta_2,beta_2s,N_bld,R_rot,*settings):

        self.Laenge = [i * 1000 for i in Element['Laenge']]
        self.DI1 = [i * 1000 for i in Element['DI1']]
        self.DI2 = [i * 1000 for i in Element['DI2']]
        self.DI3 = [i * 1000 for i in Element['DI3']]
        self.DA1 = [i * 1000 for i in Element['DA1']]
        self.DA2 = [i * 1000 for i in Element['DA2']]
        self.DA3 = [i * 1000 for i in Element['DA3']]

        self.elem_type1 = Element['elem_type1']
        self.elem_type2 = Element['elem_type2']
        self.elem_type3 = Element['elem_type3']

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

        if 'auto_rotor' in settings:
            # Automatically defining the radius of the rotor from the pickle file
            self.R_rot = (Element['parameters']['comp1']['Rrot'])*1000

        if 'manual_rotor' in settings:
            self.R_rot = R_rot

        else:
            self.R_rot = (Element['parameters']['comp1']['Rrot'])*1000
            
        #---------------------OR
            '''
            for k in range(len(self.Laenge)):
                if self.elem_type1[k]!='COMP1' or self.elem_type2[k]!='COMP1' or self.elem_type3[k]!='COMP1':
                            break
            self.R_rot = self.DA3[k]/2
            '''

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

        return

    # Defining a method to model the hub
    def hub(self,*settings):

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])
        y_5 = np.full_like(x_5, self.R_rot)  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h + self.c:
                self.coords_hub_before_3.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)
        
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3, includeCurrent=False)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0)
            .close())



        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0),-90)
                .rotate((0,0,0),(0,0,1),90)) # a Workplane object
        
        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    # Defining a method to model the hub
    def hub_v2(self,*settings):

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])
        y_5 = np.full_like(x_5, self.R_rot)  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)
        
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split1, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0)
            .close())



        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0),-90)
                .rotate((0,0,0),(0,0,1),90)) # a Workplane object
        
        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    # Defining a method to model the hub
    def hub_v3(self,*settings):

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub) #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])
        y_5 = np.full_like(x_5, self.R_rot)  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)

        tangent_12 = self.compute_tangent_between_splines(self.coords_hub_before_3_split1, self.coords_hub_before_3_split2)
        # Create tangents list for the first spline
        # All tangents are None except for the last one
        tangents_first_spline = [None] * (len(self.coords_hub_before_3_split1) - 1) + [tangent_12]

        # Tangents for the second spline
        # The first tangent is the same as the last tangent of the first spline
        tangents_second_spline = [tangent_12] + [None] * (len(self.coords_hub_before_3_split2))
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split1, tangents=tangents_first_spline, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, tangents=tangents_second_spline, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0)
            .close())



        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0),-90)
                .rotate((0,0,0),(0,0,1),90)) # a Workplane object
        
        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    def hub_v4(self,*settings):

        epsilon = 0 #mm

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])+epsilon
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub)+epsilon #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)+epsilon  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])+epsilon
        y_5 = np.full_like(x_5, self.R_rot)+epsilon  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)

        tangent_12 = self.compute_tangent_between_splines(self.coords_hub_before_3_split1, self.coords_hub_before_3_split2)
        # Create tangents list for the first spline
        # All tangents are None except for the last one
        tangents_first_spline = [None] * (len(self.coords_hub_before_3_split1) - 1) + [tangent_12]

        # Tangents for the second spline
        # The first tangent is the same as the last tangent of the first spline
        tangents_second_spline = [tangent_12] + [None] * (len(self.coords_hub_before_3_split2))
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split1, tangents=tangents_first_spline, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, tangents=tangents_second_spline, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0)
            .close())



        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0),-90)
                .rotate((0,0,0),(0,0,1),90)) # a Workplane object
        
        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    def hub_v5(self,*settings):

        epsilon = 0.01 #mm

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])+epsilon
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub)+epsilon #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)+epsilon  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])+epsilon
        y_5 = np.full_like(x_5, self.R_rot)+epsilon  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)

        tangent_12 = self.compute_tangent_between_splines(self.coords_hub_before_3_split1, self.coords_hub_before_3_split2)
        # Create tangents list for the first spline
        # All tangents are None except for the last one
        tangents_first_spline = [None] * (len(self.coords_hub_before_3_split1) - 1) + [tangent_12]

        # Tangents for the second spline
        # The first tangent is the same as the last tangent of the first spline
        tangents_second_spline = [tangent_12] + [None] * (len(self.coords_hub_before_3_split2))
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split1, tangents=tangents_first_spline, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, tangents=tangents_second_spline, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0+epsilon)
            .close())



        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(180,(0,0,0),(1,0,0))
                .mirror("XY", union=True)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0),-90)
                .rotate((0,0,0),(0,0,1),90)) # a Workplane object
        
        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    def hub_v6(self,*settings):

        epsilon = self.r_2h/2 #mm

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub) #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])
        y_5 = np.full_like(x_5, self.R_rot)  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] < self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h <= self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)


        # Create the hemisphere at the start
        hemisphere = (cq.Workplane("YZ")
                      .sphere(self.r_2h)
                      .translate((self.r_2h, 0, 0)))
        
        # Create the cylindrical section, starting from the top of the hemisphere
        cylinder_height = self.L_imp + self.R_rot/3 - self.r_2h  # Height of the cylinder
        cylinder = (cq.Workplane("YZ")
                    .circle(epsilon)  # Radius of the cylinder is epsilon
                    .extrude(cylinder_height))
        # Move the cylinder up so that it starts where the hemisphere ends
        cylinder = cylinder.translate((self.r_2h, 0, 0)).val()
        print('cylinder',cylinder)
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split2, includeCurrent=False)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0+epsilon/2)
            .hLineTo(0+self.r_2h)
            .close()).revolve(360,(0,0,0),(1,0,0)).translate((0, 0, 0)).val()
        print('self.sketch_hub ',self.sketch_hub )
        
        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]


        # Combine the hemisphere, cylinder, and complex shape
        hub = ((hemisphere.union(cylinder, clean=True, glue=False).union(self.sketch_hub, clean=True, glue=False))
               .translate((-abs(self.L_imp-self.shift),0,0))
               .rotate((0,0,0),(0,1,0),-90)
               .rotate((0,0,0),(0,0,1),90))


        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    def hub_v7(self,*settings):

        epsilon = 0.05*(self.r_2h/(2.582-0.05))

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])+epsilon
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub)+epsilon #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)+epsilon  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])#+epsilon
        y_5 = np.full_like(x_5, self.R_rot)#+epsilon  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        all_coordinates = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            all_coordinates.append((self.x_hub[i], self.y_hub[i]))
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)

        tangent_12 = self.compute_tangent_between_splines(self.coords_hub_before_3_split1, self.coords_hub_before_3_split2)
        # Create tangents list for the first spline
        # All tangents are None except for the last one
        tangents_first_spline = [None] * (len(self.coords_hub_before_3_split1) - 1) + [tangent_12]

        # Tangents for the second spline
        # The first tangent is the same as the last tangent of the first spline
        tangents_second_spline = [tangent_12] + [None] * (len(self.coords_hub_before_3_split2))
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split1, tangents=tangents_first_spline, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, tangents=tangents_second_spline, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0+epsilon)
            .close())
        
        # Gather all coordinates
        print('all_coordinates', all_coordinates)
        
        # Create the cylindrical section, starting from the top of the hemisphere
        cylinder_height = all_coordinates[-1][0] #sum([tup[0] for tup in all_coordinates]) #self.L_imp + self.R_rot/3  # Height of the cylinder
        cylinder = (cq.Workplane("YZ")
                    .circle(epsilon)  # Radius of the cylinder is epsilon
                    .extrude(cylinder_height))
        # Move the cylinder up so that it starts where the hemisphere ends
        cylinder = cylinder.translate((0, 0, 0)).val()

        # Combine hub and hemisphere
        self.sketch_hub = self.sketch_hub.union(cylinder, clean=True, glue=False)

        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(180,(0,0,0),(1,0,0))
                .mirror("XY", union=True)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0), -90)
                .rotate((0,0,0),(0,0,1), 90)) # a Workplane object

        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    def hub_v8(self,*settings):

        epsilon = 0.05*(self.r_2h/(2.582-0.05))# 1e-3 #mm

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])+epsilon
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub)+epsilon #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)+epsilon  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])+epsilon
        y_5 = np.full_like(x_5, self.R_rot)+epsilon  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)

        # Remove the first point from coords_hub_before_3_split1 and retrieve its y-value
        removed_point_x0, removed_point_y0 = self.coords_hub_before_3_split1.pop(0)
        removed_point_x1, removed_point_y1 = self.coords_hub_before_3_split1[0]
        print('removed_point_x0',removed_point_x0)
        print('removed_point_x1',removed_point_x1)
        print('removed_point_y0',removed_point_y0)
        print('removed_point_y1',removed_point_y1)

        tangent_12 = self.compute_tangent_between_splines(self.coords_hub_before_3_split1, self.coords_hub_before_3_split2)
        # Create tangents list for the first spline
        # All tangents are None except for the last one
        tangents_first_spline = [None] * (len(self.coords_hub_before_3_split1) - 1) + [tangent_12]

        # Tangents for the second spline
        # The first tangent is the same as the last tangent of the first spline
        tangents_second_spline = [tangent_12] + [None] * (len(self.coords_hub_before_3_split2))
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split1, tangents=tangents_first_spline, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, tangents=tangents_second_spline, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0+removed_point_y0+removed_point_y1)
            .close())
        
        # Create the cylindrical section, starting from the top of the hemisphere
        cylinder_eps_r = 0 #mm
        cylinder_height = self.L_imp + self.R_rot/3 - removed_point_x0 - removed_point_x1 # Height of the cylinder
        cylinder = (cq.Workplane("YZ")
                    .circle(removed_point_y0+removed_point_y1+cylinder_eps_r)  # Radius of the cylinder is epsilon
                    .extrude(cylinder_height))
        # Move the cylinder up so that it starts where the hemisphere ends
        cylinder = cylinder.translate((removed_point_x0 + removed_point_x1, 0, 0)).val()

        # Combine hub and hemisphere
        self.sketch_hub = self.sketch_hub.union(cylinder, clean=True, glue=False)

        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                #.mirror("XY", union=True)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0), -90)
                .rotate((0,0,0),(0,0,1), 90)) # a Workplane object

        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    def hub_v9(self,*settings):

        epsilon = 0# 1e-3 #mm
        epsilon_tip = 1#mm

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])+epsilon
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub)+epsilon #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)+epsilon+epsilon_tip  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])+epsilon
        y_5 = np.full_like(x_5, self.R_rot)+epsilon  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)

        # Remove the first point from coords_hub_before_3_split1 and retrieve its y-value
        #removed_point_x0, removed_point_y0 = self.coords_hub_before_3_split1.pop(0)
        #removed_point_x1, removed_point_y1 = self.coords_hub_before_3_split1[0]
        #print('removed_point_x0',removed_point_x0)
        ##print('removed_point_x1',removed_point_x1)
        #print('removed_point_y0',removed_point_y0)
        #print('removed_point_y1',removed_point_y1)

        tangent_12 = self.compute_tangent_between_splines(self.coords_hub_before_3_split1, self.coords_hub_before_3_split2)
        # Create tangents list for the first spline
        # All tangents are None except for the last one
        tangents_first_spline = [None] * (len(self.coords_hub_before_3_split1) - 1) + [tangent_12]

        # Tangents for the second spline
        # The first tangent is the same as the last tangent of the first spline
        tangents_second_spline = [tangent_12] + [None] * (len(self.coords_hub_before_3_split2))
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3_split1, tangents=tangents_first_spline, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, tangents=tangents_second_spline, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0)#+removed_point_y0+removed_point_y1)
            .close())

        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                #.mirror("XY", union=True)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0), -90)
                .rotate((0,0,0),(0,0,1), 90)) # a Workplane object

        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub
    
    def hub_v10(self,*settings):

        epsilon = 0.1# 1e-3 #mm

        if 'section view' in settings:
            bottom_hub = False
        else:
            bottom_hub = True


        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        n_z = 200
        # Angle for ellipse description
        alpha_hub = np.linspace(3 * np.pi / 2, 2 * np.pi, n_z, endpoint=False)
        # Define z_hub
        z_hub = lambda theta: self.r_2h + self.c * np.cos(theta)
        # Define r_hub
        r_hub = lambda theta: self.r_4 + self.d * np.sin(theta)

        # Define the intervals
        n_points_1 = round(self.r_2h / self.delta_x)
        n_points_2 = round(self.c / self.delta_x)
        n_points_3 = round(self.b_6 / self.delta_x)
        n_points_4 = round((self.L_imp - self.r_2h - self.c - self.b_6) / self.delta_x)
        n_points_5 = round((self.R_rot/3) / self.delta_x)

        x_1 = np.linspace(0, self.r_2h, n_points_1, endpoint=False)
        x_2 = z_hub(alpha_hub) #np.linspace(self.r_2h, self.r_2h + self.c, n_points_2, endpoint=False)
        x_3 = np.linspace(self.r_2h + self.c, self.r_2h + self.c + self.b_6, n_points_3, endpoint=False)
        x_4 = np.linspace(self.r_2h + self.c + self.b_6, self.L_imp, n_points_4, endpoint=False)
        x_5 = np.linspace(self.L_imp, self.L_imp + self.R_rot/3, n_points_5, endpoint=True)

        # Concatenate all x arrays
        self.x_hub = np.concatenate((x_1, x_2, x_3, x_4, x_5))

        # Apply the piecewise function
        ##y_1 = np.sqrt(self.r_2h**2 - (x_1 - self.r_2h)**2)
        y_1 = np.array([np.sqrt(self.r_2h**2 - (x_val - self.r_2h)**2) if (self.r_2h**2 - (x_val - self.r_2h)**2) >= 0 else 0 for x_val in x_1])+epsilon
        ##y_2 = self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_2 - self.r_2h)**2) 
        y_2 = r_hub(alpha_hub)+epsilon #np.array([self.r_4 - np.sqrt(self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) if (self.d**2 - (self.d**2/self.c**2)*(x_val - self.r_2h)**2) >= 0 else self.r_4 for x_val in x_2])
        y_3 = np.full_like(x_3, self.r_4)+epsilon  # create an array filled with r_4 value
        ##y_4 = self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_4 - self.L_imp)**2)
        y_4 = np.array([self.r_4 - np.sqrt(self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) if (self.f**2 - (self.f**2/self.e**2)*(x_val - self.L_imp)**2) >= 0 else self.r_4 for x_val in x_4])+epsilon
        y_5 = np.full_like(x_5, self.R_rot)+epsilon  # create an array filled with R_rot value


        # Concatenate all y arrays
        self.y_hub = np.concatenate((y_1, y_2, y_3, y_4, y_5))

        print('self.x_hub',self.x_hub)
        print('self.y_hub',self.y_hub)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3_split1 = []
        self.coords_hub_before_3_split2 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Create tuple coordinates and segregate to respective parts
        for i in range(len(self.x_hub)):
            if self.x_hub[i] <= self.r_2h:
                self.coords_hub_before_3_split1.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h < self.x_hub[i] <= self.r_2h + self.c :
                self.coords_hub_before_3_split2.append((self.x_hub[i], self.y_hub[i]))
            elif self.r_2h + self.c < self.x_hub[i] <= self.r_2h + self.c + self.b_6:
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
            else:
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))

        print('self.coords_hub_before_3',self.coords_hub_before_3_split1)
        print('self.coords_hub_before_3',self.coords_hub_before_3_split2)
        print('self.coords_hub_3',self.coords_hub_3)
        print('self.coords_hub_after_3',self.coords_hub_after_3)

        # Remove the first point from coords_hub_before_3_split1 and retrieve its y-value
        removed_point_x0, removed_point_y0 = self.coords_hub_before_3_split1.pop(0)
        removed_point_x1, removed_point_y1 = self.coords_hub_before_3_split1[0]
        print('removed_point_x0',removed_point_x0)
        print('removed_point_x1',removed_point_x1)
        print('removed_point_y0',removed_point_y0)
        print('removed_point_y1',removed_point_y1)

        tangent_12 = self.compute_tangent_between_splines(self.coords_hub_before_3_split1, self.coords_hub_before_3_split2)
        # Create tangents list for the first spline
        # All tangents are None except for the last one
        tangents_first_spline = [None] * (len(self.coords_hub_before_3_split1) - 1) + [tangent_12]

        # Tangents for the second spline
        # The first tangent is the same as the last tangent of the first spline
        tangents_second_spline = [tangent_12] + [None] * (len(self.coords_hub_before_3_split2))
        
        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .vLineTo(0+epsilon)
            .spline(self.coords_hub_before_3_split1, tangents=tangents_first_spline, includeCurrent=False)
            .spline(self.coords_hub_before_3_split2, tangents=tangents_second_spline, includeCurrent=True)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0)
            .close())
        

        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        # NO SPLIT
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                #.mirror("XY", union=True)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0), -90)
                .rotate((0,0,0),(0,0,1), 90)) # a Workplane object

        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turbocompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly, hub

    def compute_tangent_between_splines(self, coords1, coords2):
        # Ensure that both lists have at least one element
        if not coords1 or not coords2:
            raise ValueError("Coordinate lists must not be empty")

        # Last point of the first spline and first point of the second spline
        last_point_coords1 = coords1[-1]
        first_point_coords2 = coords2[0]

        # Compute the tangent vector
        tangent_vector = (first_point_coords2[0] - last_point_coords1[0], 
                        first_point_coords2[1] - last_point_coords1[1])

        # Normalizing the tangent vector (optional, depends on your use case)
        magnitude = (tangent_vector[0]**2 + tangent_vector[1]**2)**0.5
        normalized_tangent = (tangent_vector[0] / magnitude, tangent_vector[1] / magnitude)

        return normalized_tangent
    
    # Defining a method to retrieve the coordinates of the blades from the excel file
    def blades_excel(self,filename):

        # Storing lists of non-planar blade curves
        points_blade = {}
        coords_blade = {}
        
        # Reading the excel file with blade points
        excelfile = self.cwf  + '/COMP/' + filename

        # Detects the type of spreadsheet inputted to use the correct reading command
        split_string = excelfile.split('.')
        last_element = split_string[-1]
        
        if last_element == 'xls':
            df_blade = pd.read_excel(io=excelfile,header=None)
        elif last_element == 'xlsx':
            df_blade = pd.read_excel(excelfile,engine='openpyxl',header=None)
        else:
            raise TypeError('Impeller.blades_excel: Provide a valid excel file.')
        
        # Setting index markers to indicate start and end of storing values
        start = None
        end = None
        
        # Looping to iterate over all the rows
        for index, row in df_blade.iterrows():
        
        # Updating index markers every time 'StartCurve' or 'EndCurve' is detected
            if 'StartCurve' in row.values:
                start = index+ 1 #skip the header row
            elif 'EndCurve' in row.values:
                end = index
        
        # Storing curve points
                if len(points_blade) == 0:
                    points_blade[0] = df_blade.iloc[start:end,:]
                else:
                    points_blade[len(points_blade)] = df_blade.iloc[start:end,:]
                start = None
                end = None

        # Grouping the curve x,y,z coordinates in tuples
        for i in range(len(points_blade)):
            coords_blade[i] = list(zip(points_blade[i].iloc[:,0].astype(float),\
                                            points_blade[i].iloc[:,1].astype(float),\
                                            points_blade[i].iloc[:,2].astype(float)))
        
        return coords_blade
    
    # Defining a method to model the blades
    def model_blades(self,coords):
        # To shift the blades further into hub and ensure connection
        epsilon = 0 #10e-6*1000 # everything is in mm basis here
        
        # Initializing lists to store the solids
        blade_solid={}

        # Sketching the first blade curve as an edge
        blade_edge0 = cq.Edge.makeSpline([cq.Vector(p) for p in coords[0]][0:-1]).close()
        
        # Creating a closed surface within the first blade curve edge
        blade_face0 = cq.Face.makeNSidedSurface([blade_edge0],[])
        
        # Sketching the last blade curve as an edge
        blade_edge_last = cq.Edge.makeSpline([cq.Vector(p) for p in coords[len(coords)-1]][0:-1]).close()
        
        # Creating a closed surface within the last blade curve edge
        blade_face_last = cq.Face.makeNSidedSurface([blade_edge_last],[])
        
        # Lofting the closed surfaces to produce the overall blade shape
        blade_lofted = cq.Solid.makeLoft([blade_edge0,blade_edge_last],True)
        
        # Shelling the faces and edges
        blade_shell = cq.Shell.makeShell([blade_face0,blade_face_last,blade_lofted]).fix()
        
        # Solidifying the produced shell and rotating
        # blade_solid[0] = cq.Solid.makeSolid(blade_shell).translate((0,0,-14))
        blade_solid[0] = cq.Solid.makeSolid(blade_shell).translate((0,0,-abs(self.L_imp-self.shift)+epsilon))
        #print('abs(self.L_imp-self.shift)',abs(self.L_imp-self.shift))

        print('type(blade_solid[0])',type(blade_solid[0]))

        return blade_solid
    
    def model_blades_v2(self, coords):
        epsilon = 0  # Adjust as needed for precision

        def calculate_tangents(points):
            tangents = []
            for i in range(len(points)):
                if i == 0:
                    # Tangent for the first point
                    tangent = points[i + 1] - points[i]
                elif i == len(points) - 1:
                    # Tangent for the last point
                    tangent = points[i] - points[i - 1]
                else:
                    # Tangent for intermediate points
                    tangent = points[i + 1] - points[i - 1]

                # Normalize the tangent vector
                tangent_length = tangent.Length
                if tangent_length != 0:
                    tangent = tangent.multiply(1.0 / tangent_length)

                tangents.append(tangent)
            return tangents

        # Function to create a spline with automatic tangent calculation
        def create_spline_with_tangents(points):
            vectors = [cq.Vector(p) for p in points]
            tangents = calculate_tangents(vectors)
            return cq.Edge.makeSpline(vectors, tangents=tangents)

        # Creating the first and last blade edges as splines
        blade_edge0 = create_spline_with_tangents(coords[0]).close()
        blade_edge_last = create_spline_with_tangents(coords[-1]).close()

        # Creating a surface within the first blade curve
        blade_face0 = cq.Face.makeNSidedSurface([blade_edge0], [])

        # Creating a surface within the last blade curve
        blade_face_last = cq.Face.makeNSidedSurface([blade_edge_last], [])

        # Lofting between the two faces
        blade_lofted = cq.Solid.makeLoft([blade_face0, blade_face_last])

        # Final transformation and solidification
        blade_solid = blade_lofted.translate((0, 0, -abs(self.L_imp - self.shift) + epsilon))

        return blade_solid
    
    

    def model_blades_v3(self, coords):
        blade_solid = {}

        # Extract the hub and shroud curves
        hub_curve, shroud_curve = coords

        # Check that the hub and shroud have the same number of points
        if len(hub_curve) != len(shroud_curve):
            raise ValueError("Hub and shroud curves must have the same number of points")

        # Create wires for each cross-section by pairing corresponding points
        section_wires = []
        num_points = len(hub_curve)
        for i in range(num_points // 2):  # We assume an even number of points
            # Pair the corresponding points from the hub and shroud
            # The first and last point of the hub pair with the first and last of the shroud, and so on.
            hub_points = [hub_curve[i], hub_curve[num_points - 1 - i]]
            shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            # Combine the points to form the cross-section wire
            cross_section_points = hub_points + shroud_points[::-1]  # Reverse shroud points for correct ordering
            section_wire = cq.Wire.makePolygon([cq.Vector(*p) for p in cross_section_points]).close()
            section_wires.append(section_wire)

        # Loft through the section wires to create the blade
        blade_lofted = cq.Solid.makeLoft(section_wires, True)

        # Translate the blade solid to the correct position
        blade_solid[0] = blade_lofted.translate((0, 0, -abs(self.L_imp - self.shift)))

        return blade_solid
    
    def model_blades_v4(self, coords):
        epsilon = 0  # Adjust as needed for precision

        blade_solid = {}

        # Extract the hub and shroud curves
        hub_curve, shroud_curve = coords

        # Create wires for each cross-section
        section_wires = []
        num_points = len(hub_curve)
        for i in range(num_points // 2):  # We assume an even number of points
            # Pair the corresponding points from the hub and shroud
            # The first and last point of the hub pair with the first and last of the shroud, and so on.
            hub_points = [hub_curve[i], hub_curve[num_points - 1 - i]]
            shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            if i == 0:
                first_hub_points    = [hub_curve[i], hub_curve[num_points - 1 - i]]
                first_shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            if i == num_points // 2 -1: 
                last_hub_points     = [hub_curve[i], hub_curve[num_points - 1 - i]]
                last_shroud_points  = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            # Combine the points to form the cross-section wire
            cross_section_points = hub_points + shroud_points[::-1]  # Reverse shroud points for correct ordering
            #print('cross_section_points',cross_section_points)
            section_wire = cq.Wire.makePolygon([cq.Vector(*p) for p in cross_section_points]).close()
            section_wires.append(section_wire)

        first_hub_curve    = cq.Edge.makeSpline([cq.Vector(p) for p in first_hub_points])
        first_shroud_curve = cq.Edge.makeSpline([cq.Vector(p) for p in first_shroud_points])

        # Creating a closed surface within the first blade curve edge
        first_face = cq.Face.makeRuledSurface(first_hub_curve, first_shroud_curve)

        last_hub_curve    = cq.Edge.makeSpline([cq.Vector(p) for p in last_hub_points])
        last_shroud_curve = cq.Edge.makeSpline([cq.Vector(p) for p in last_shroud_points])

        # Creating a closed surface within the first blade curve edge
        last_face  = cq.Face.makeRuledSurface(last_hub_curve, last_shroud_curve)

        # Loft through the section wires to create the blade surface
        blade_lofted = cq.Solid.makeLoft(section_wires, True)

        # To ensure a solid, create a shell from the lofted surfaces
        blade_shell = cq.Shell.makeShell([first_face, blade_lofted,last_face]).fix()

        # Solidify the shell to create a solid blade
        blade_solidified = cq.Solid.makeSolid(blade_shell)

        # Translate the blade solid to the correct position
        blade_solid[0] = blade_solidified.translate((0, 0, -abs(self.L_imp - self.shift) + epsilon))

        return blade_solid
    
    def model_blades_v5(self, coords):
        epsilon = 0  # Adjust as needed for precision

        blade_solid = {}

        # Extract the hub and shroud curves
        hub_curve, shroud_curve = coords

        # Create wires for each cross-section
        section_wires = []
        num_points = len(hub_curve)
        for i in range(num_points // 2):  # We assume an even number of points
            # Pair the corresponding points from the hub and shroud
            hub_points = [hub_curve[i], hub_curve[num_points - 1 - i]]
            shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            
            # Combine the points to form the cross-section wire
            cross_section_points = hub_points + shroud_points[::-1]  # Reverse shroud points for correct ordering
            section_wire = cq.Wire.makePolygon([cq.Vector(*p) for p in cross_section_points]).close()
            section_wires.append(section_wire)

        # Loft through the section wires to create the blade surface
        blade_lofted = cq.Solid.makeLoft(section_wires, ruled=True).Faces()

        # Sketching the first blade curve as an edge
        blade_edge0 = cq.Edge.makeSpline([cq.Vector(p) for p in coords[0]][0:-1]).close()
        
        # Creating a closed surface within the first blade curve edge
        end_face1 = cq.Face.makeNSidedSurface([blade_edge0],[])
        
        # Sketching the last blade curve as an edge
        blade_edge_last = cq.Edge.makeSpline([cq.Vector(p) for p in coords[len(coords)-1]][0:-1]).close()
        
        # Creating a closed surface within the last blade curve edge
        end_face2 = cq.Face.makeNSidedSurface([blade_edge_last],[])

        # Create planar end faces from the first and last wires
        #end_face1 = cq.Face.makeNSidedSurface(section_wires[0],[])
        #end_face2 = cq.Face.makeNSidedSurface(section_wires[-1],[])

        # Create a shell by combining the loft with the end faces
        blade_shell = cq.Shell.makeShell(blade_lofted+[end_face1,end_face2]).fix()

        # Solidify the shell to create a solid blade
        blade_solidified = cq.Solid.makeSolid(blade_shell)

        # Translate the blade solid to the correct position
        blade_solid[0] = blade_solidified.translate((0, 0, -abs(self.L_imp - self.shift) + epsilon))

        return blade_solid
    
    def model_blades_v6(self, coords):
        epsilon = 0  # Adjust as needed for precision

        blade_solid = {}

        # Extract the hub and shroud curves
        hub_curve, shroud_curve = coords

        # Create wires for each cross-section
        section_wires     = []
        section_surface   = []
        for i in range(len(hub_curve)):
            points = [hub_curve[i], shroud_curve[i]]  # Pair points to form edges
            # More logic to ensure correct point pairing goes here
            wire = cq.Wire.makePolygon([cq.Vector(*p) for p in points]).close()
            section_wires.append(wire)

            surface = cq.Face.makeNSidedSurface([wire], [])
            section_surface.append(surface)


        # Perform the loft operation
        blade_lofted = cq.Solid.makeLoft(section_wires)

        # Create end wires from the first and last cross-sections
        end_faces = [section_surface[0],section_surface[-1]]

        # Combine the lofted blade with the end faces to form a shell
        blade_faces = [blade_lofted] + end_faces
        blade_shell = cq.Shell.makeShell(blade_faces).fix()

        # Solidify the shell to create a solid blade
        blade_solidified = cq.Solid.makeSolid(blade_shell)

        # Translate the blade solid to the correct position
        blade_solid[0] = blade_solidified.translate((0, 0, -abs(self.L_imp - self.shift) + epsilon))

        return blade_solid
    
    def model_blades_v7(self, coords):
        blade_solid = {}

        epsilon = self.b_6#*(2/3)

        # Extract the hub and shroud curves
        hub_curve, shroud_curve = coords

        # Create wires for each cross-section
        section_wires = []
        num_points = len(hub_curve)

        # Extend hub points into hub:
        for i in range(num_points):
            # Calculate the vector from the shroud point to the hub point
            direction_vector = cq.Vector(hub_curve[i]) - cq.Vector(shroud_curve[i])
            # Scale the vector by epsilon to get the extension vector
            extension_vector = direction_vector.normalized() * epsilon
            # Extend the hub point
            new_hub_point = cq.Vector(hub_curve[i]) + extension_vector
            hub_curve[i] = (new_hub_point.x, new_hub_point.y, new_hub_point.z)

        for i in range(num_points // 2):

            hub_points = [hub_curve[i], hub_curve[num_points - 1 - i]]
            shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            if i == 0:
                first_hub_points    = [hub_curve[i], hub_curve[num_points - 1 - i]]
                first_shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            if i == num_points // 2 -1: 
                last_hub_points     = [hub_curve[i], hub_curve[num_points - 1 - i]]
                last_shroud_points  = [shroud_curve[i], shroud_curve[num_points - 1 - i]]

            # Combine the points to form the cross-section wire
            cross_section_points = hub_points + shroud_points[::-1] # Reverse shroud points for correct ordering
            section_wire = cq.Wire.makePolygon([cq.Vector(p) for p in cross_section_points]).close()
            section_wires.append(section_wire)


        first_hub_curve    = cq.Edge.makeSpline([cq.Vector(p) for p in first_hub_points])
        first_shroud_curve = cq.Edge.makeSpline([cq.Vector(p) for p in first_shroud_points])

        # Creating a closed surface within the first blade curve edge
        first_face = cq.Face.makeRuledSurface(first_hub_curve, first_shroud_curve)

        last_hub_curve    = cq.Edge.makeSpline([cq.Vector(p) for p in last_hub_points])
        last_shroud_curve = cq.Edge.makeSpline([cq.Vector(p) for p in last_shroud_points])

        # Creating a closed surface within the first blade curve edge
        last_face  = cq.Face.makeRuledSurface(last_hub_curve, last_shroud_curve)

        # Loft through the section wires to create the blade surface
        blade_lofted = cq.Solid.makeLoft(section_wires, True)

        # To ensure a solid, create a shell from the lofted surfaces
        blade_shell = cq.Shell.makeShell([first_face, blade_lofted,last_face]).fix()

        # Solidify the shell to create a solid blade
        blade_solidified = cq.Solid.makeSolid(blade_shell)

        # Translate the blade solid to the correct position
        blade_solid[0] = blade_solidified.translate((0, 0, -abs(self.L_imp - self.shift)))

        return blade_solid
    
    def model_blades_v8(self, coords):
        blade_solid = {}

        epsilon = self.b_6#*(2/3)

        # Extract the hub and shroud curves
        hub_curve, shroud_curve = coords

        # Initialize variables to store the spline points
        curve_hub_1_points = []
        curve_hub_2_points = []
        curve_shroud_1_points = []
        curve_shroud_2_points = []

        # Create wires for each cross-section
        section_wires = []
        num_points = len(hub_curve)

        # Extend hub points into hub:
        for i in range(num_points):
            # Calculate the vector from the shroud point to the hub point
            direction_vector = cq.Vector(hub_curve[i]) - cq.Vector(shroud_curve[i])
            # Scale the vector by epsilon to get the extension vector
            extension_vector = direction_vector.normalized() * epsilon
            # Extend the hub point
            new_hub_point = cq.Vector(hub_curve[i]) + extension_vector
            hub_curve[i] = (new_hub_point.x, new_hub_point.y, new_hub_point.z)

        for i in range(num_points // 2):

            hub_points = [hub_curve[i], hub_curve[num_points - 1 - i]]
            shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            if i == 0:
                first_hub_points    = [hub_curve[i], hub_curve[num_points - 1 - i]]
                first_shroud_points = [shroud_curve[i], shroud_curve[num_points - 1 - i]]
            if i == num_points // 2 -1: 
                last_hub_points     = [hub_curve[i], hub_curve[num_points - 1 - i]]
                last_shroud_points  = [shroud_curve[i], shroud_curve[num_points - 1 - i]]

            curve_hub_1_points.append(hub_curve[i])
            curve_hub_2_points.append(hub_curve[num_points - 1 - i])

            curve_shroud_1_points.append(shroud_curve[i])
            curve_shroud_2_points.append(shroud_curve[num_points - 1 - i])

            # Combine the points to form the cross-section wire
            cross_section_points = hub_points + shroud_points[::-1] # Reverse shroud points for correct ordering
            section_wire = cq.Wire.makePolygon([cq.Vector(p) for p in cross_section_points]).close()
            section_wires.append(section_wire)

        # Creating the splines to build the 6 rules faces of the blade
        first_hub_curve    = cq.Edge.makeSpline([cq.Vector(p) for p in first_hub_points])
        first_shroud_curve = cq.Edge.makeSpline([cq.Vector(p) for p in first_shroud_points])

        last_hub_curve    = cq.Edge.makeSpline([cq.Vector(p) for p in last_hub_points])
        last_shroud_curve = cq.Edge.makeSpline([cq.Vector(p) for p in last_shroud_points])

        curve_hub_1    = cq.Edge.makeSpline([cq.Vector(p) for p in curve_hub_1_points])
        curve_hub_2    = cq.Edge.makeSpline([cq.Vector(p) for p in curve_hub_2_points])
        curve_shroud_1 = cq.Edge.makeSpline([cq.Vector(p) for p in curve_shroud_1_points])
        curve_shroud_2 = cq.Edge.makeSpline([cq.Vector(p) for p in curve_shroud_2_points])

        # Creating a closed surface within the first blade curve edge (blade face 0)
        first_face = cq.Face.makeRuledSurface(first_hub_curve, first_shroud_curve)

        # Creating a closed surface within the first blade curve edge (blade face 5)
        last_face  = cq.Face.makeRuledSurface(last_hub_curve, last_shroud_curve)

        # Creating remaining closed surfaces of the blade
        second_face  = cq.Face.makeRuledSurface(curve_hub_1, curve_shroud_1)
        third_face   = cq.Face.makeRuledSurface(curve_shroud_1, curve_shroud_2)
        fourth_face  = cq.Face.makeRuledSurface(curve_shroud_2, curve_hub_2)
        fifth_face   = cq.Face.makeRuledSurface(curve_hub_2, curve_hub_1)

        # List of faces
        list_faces = [first_face, second_face, third_face, fourth_face, fifth_face, last_face] 

        # To ensure a solid, create a shell from the lofted surfaces
        blade_shell = cq.Shell.makeShell(list_faces).fix()

        # Solidify the shell to create a solid blade
        blade_solidified = cq.Solid.makeSolid(blade_shell)

        # Translate the blade solid to the correct position
        blade_solid[0] = blade_solidified.translate((0, 0, -abs(self.L_imp - self.shift)))

        return blade_solid
    
    def trim_blades(self, blade):

        # Create a larger outer cylinder
        outer_cylinder = cq.Solid.makeCylinder(2*self.r_4, self.L_imp, cq.Vector(0, 0, -abs(self.L_imp-self.shift)))

        # Create a smaller inner cylinder
        inner_cylinder = cq.Solid.makeCylinder(self.r_4, self.L_imp, cq.Vector(0, 0, -abs(self.L_imp-self.shift)))

        # Subtract the smaller cylinder from the larger one to create a hollow cylinder
        hollow_cylinder = outer_cylinder.cut(inner_cylinder, clean=False)

        #print('blade[0]',blade[0])

        # Use a boolean cut operation to subtract the hollow cylinder from the blade
        ##blade[0] = blade[0].cut(hollow_cylinder)

        #print('blade[0]',blade[0])

        # Create a Workplane
        ##workplane = cq.Workplane("XY")

        # Add the blade to the Workplane
        ##workplane = workplane.add(blade[0])

        # Perform the union operation
        ##result = workplane.intersect(hollow_cylinder)

        # Update the blade[0] with the result
        ##blade[0] = result.val()
        blade_tmp = blade[0]
        blade[0] = blade_tmp.cut(hollow_cylinder) #.copy()
        

        return blade
    
    def trim_blades_old(self, blade):

        # Create a Workplane
        workplane = cq.Workplane("XY")

        # Add the blade to the Workplane
        workplane = workplane.add(blade[0])

        # Create a larger outer cylinder
        outer_cylinder = cq.Workplane("XY").circle(2*self.r_4).extrude(self.L_imp)

        # Create a smaller inner cylinder
        inner_cylinder = cq.Workplane("XY").circle(self.r_4).extrude(self.L_imp)

        # Subtract the smaller cylinder from the larger one to create a hollow cylinder
        hollow_cylinder = outer_cylinder.cut(inner_cylinder)

        # Perform the union operation
        result = workplane.cut(hollow_cylinder)

        # Update the blade[0] with the result
        blade[0] = result.val()

        

        return blade

    '''
    # Defining a method to pattern the blades
    def rotate_blade(self,blade,bladename):

        assembly = cq.Assembly(name = bladename)
        assembly.add(blade[0],color=cq.Color('red3'), name = bladename + ' 1')

        # Looping to model the blade pattern
        for i in range(0,self.N_bld-1):
        
            # Rotating about the x axis by the corresponding angle
            blade[i+1] = (blade[i].transformed((0,0,360/self.N_bld)))
            assembly.add(blade[i+1],color=cq.Color('red3'), name = bladename + ' ' + str(i+2))

        return assembly
    '''
    # Defining a method to pattern the blades
    def rotate_blade(self, blade, bladename):

        blades_solid = []
        blades_solid.append(blade[0])

        assembly = cq.Assembly(name=bladename)
        assembly.add(blade[0], color=cq.Color('red3'), name=bladename + ' 1')

        # Looping to model the blade pattern
        for i in range(0, self.N_bld-1):
        
            # Rotating about the Z axis by the corresponding angle
            # The rotation center is assumed to be (0, 0, 0) for simplicity
            blade[i+1] = blade[i].rotate((0, 0, 0), (0, 0, 1), 360/self.N_bld)
            blades_solid.append(blade[i+1])
            assembly.add(blade[i+1], color=cq.Color('red3'), name=bladename + ' ' + str(i+2))



        return assembly, blades_solid
    
    
    # Defining a method to combine the impeller components in a common assembly
    def assemble(self,files,*settings):
        assembly = cq.Assembly(name='Compressor')

        for i in range(0,len(files)):
            assembly.add(files[i],name='Subassembly '+str(i+1))

        #if 'stl' or 'STL' in settings:
        #    cq.exporters.export(assembly.toCompound(), self.cwf + '/STL/Compressor.stl')

        assembly.save(self.cwf  + '/STEP/Compressor.step')

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    def assemble_v2(self, components, *settings):
        # components is expected to be a tuple containing (Hub, Mainblades, Splitterblades)
        # Extracting hub and blade solids from the assemblies
        hub_solid = components[0]  # Assuming Hub is a solid
        main_blades_assembly = components[1]
        splitter_blades_assembly = components[2]

        # Extract solids from blade assemblies
        main_blade_solids = [item for item in main_blades_assembly.children] #if isinstance(item, cq.Solid)]
        print('main_blade_solids',main_blade_solids)
        splitter_blade_solids = [item for item in splitter_blades_assembly.children] #if isinstance(item, cq.Solid)]
        print('splitter_blade_solids',splitter_blade_solids)

        # Union the main blades with the hub
        for blade in main_blade_solids:
            hub_solid = hub_solid.union(blade)

        # Union the splitter blades with the hub
        for blade in splitter_blade_solids:
            hub_solid = hub_solid.union(blade)


        # Rest of your assembly code...
        # For example, save the combined object as a STEP or STL file
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='Combined Impeller')
        # ...

        assembly.save(self.cwf  + '/STEP/Compressor.step')

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    def assemble_v3(self, components, *settings):
       # Extracting the hub and blade solids from the components
        hub_solid = components[0]  # The hub as a solid
        blade_solids = components[1]  # List of blade solids
        splitter_blade_solids = components[2]  # List of splitter blade solids
        print('hub_solid',hub_solid)
        print('blade_solids',blade_solids)
        print('splitter_blade_solids',splitter_blade_solids)
        # Union the blades with the hub
        for blade in blade_solids:
            hub_solid = hub_solid.union(blade, clean=True, glue=False)
            print('passed a blade union with hub')
        
        # Union the splitter blades with the hub
        for splitter_blade in splitter_blade_solids:
            hub_solid = hub_solid.union(splitter_blade,clean=True, glue=False)
            print('passed a splitter union with hub')
        

        # Rest of the assembly code, like saving the combined object
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='Combined Impeller')

        assembly.save(self.cwf  + '/STEP/Compressor.step')
        print('Pausing 5 seconds for writing Comrpessor STEP.')
        time.sleep(5)

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly

    def assemble_v4(self, components, *settings):
       # Extracting the hub and blade solids from the components
        hub_solid = components[0]  # The hub as a solid
        blade_solids = components[1]  # List of blade solids
        splitter_blade_solids = components[2]  # List of splitter blade solids
        print('hub_solid',hub_solid)
        print('blade_solids',blade_solids)
        print('splitter_blade_solids',splitter_blade_solids)
        # Union the blades with the hub
        for blade in blade_solids:
            hub_solid = hub_solid+blade
            print('passed a blade union with hub')
        
        # Union the splitter blades with the hub
        for splitter_blade in splitter_blade_solids:
            hub_solid = hub_solid+splitter_blade
            print('passed a splitter union with hub')
        

        # Rest of the assembly code, like saving the combined object
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='Combined Impeller')

        assembly.save(self.cwf  + '/STEP/Compressor.step')

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    def assemble_v5(self, components, *settings):
       # Extracting the hub and blade solids from the components
        hub_solid = components[0].val()#.Solids()[0]  # The hub as a solid
        blade_solids = components[1]  # List of blade solids
        splitter_blade_solids = components[2]  # List of splitter blade solids
        print('hub_solid',hub_solid)
        print('blade_solids',blade_solids)
        print('splitter_blade_solids',splitter_blade_solids)
        # Union the blades with the hub
        for blade in blade_solids:
            #hub_solid = hub_solid.union(blade, clean=True, glue=False)
            hub_solid = hub_solid.fuse(blade)
            print('passed a blade fuse with hub')
        
        # Union the splitter blades with the hub
        for splitter_blade in splitter_blade_solids:
            #hub_solid = hub_solid.union(splitter_blade,clean=True, glue=False)
            hub_solid = hub_solid.fuse(splitter_blade)
            print('passed a splitter fuse with hub')
        

        # Rest of the assembly code, like saving the combined object
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='Combined Impeller')

        assembly.save(self.cwf  + '/STEP/Compressor.step')

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    def assemble_v6(self, components, *settings):
       # Extracting the hub and blade solids from the components
        hub_solid = components[0].val()#.Solids()[0]  # The hub as a solid
        blade_solids = components[1]  # List of blade solids
        splitter_blade_solids = components[2]  # List of splitter blade solids
        print('hub_solid',hub_solid)
        print('blade_solids',blade_solids)
        print('splitter_blade_solids',splitter_blade_solids)

        cq.exporters.export(hub_solid,self.cwf  + '/STEP/hub_solid.step', opt={"write_pcurves": False, "precision_mode": -1})
        cq.exporters.export(blade_solids[0],self.cwf  + '/STEP/blade_solids.step', opt={"write_pcurves": False, "precision_mode": -1})
        cq.exporters.export(splitter_blade_solids[0],self.cwf  + '/STEP/splitter_blade_solids.step', opt={"write_pcurves": False, "precision_mode": -1})
        
        hub_solid= cq.importers.importStep(self.cwf  + '/STEP/hub_solid.step')
        blade_solids = cq.importers.importStep(self.cwf  + '/STEP/blade_solids.step')
        splitter_blade_solids = cq.importers.importStep(self.cwf  + '/STEP/splitter_blade_solids.step')

        # Rest of the assembly code, like saving the combined object
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='hub_solid')
        assembly.add(blade_solids, name='blade_solids')
        assembly.add(splitter_blade_solids, name='splitter_blade_solids')

        assembly.save(self.cwf  + '/STEP/Compressor.step', mode="fused", clean = True, glue=False, write_pcurves=False)
        print('Pausing 5 seconds for writing Comrpessor STEP.')
        time.sleep(5)

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    def assemble_v7(self, components, *settings):
       # Extracting the hub and blade solids from the components
        hub_solid = components[0]  # The hub as a solid
        blade_solids = components[1]  # List of blade solids
        splitter_blade_solids = components[2]  # List of splitter blade solids
        print('hub_solid',hub_solid)
        print('blade_solids',blade_solids)
        print('splitter_blade_solids',splitter_blade_solids)
        # Union the blades with the hub
        self.N_bld

        for i in range(self.N_bld):
            if i == 0:
                blade = blade_solids[i]
            else:
                blade = blade.fuse(blade_solids[i])
                print('passed a blade fused')
        
        # Union the splitter blades with the hub
        for i in range(self.N_bld):
            if i == 0:
                splitter_blade = splitter_blade_solids[i]
            else:
                splitter_blade = splitter_blade.fuse(splitter_blade_solids[i])
                print('passed a splitter fused')
            
        blades = blade.fuse(splitter_blade)
        hub_solid = hub_solid.union(blades,clean=True, glue=False)
        print('passed a union with hub')

        # Rest of the assembly code, like saving the combined object
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='Combined Impeller')

        assembly.save(self.cwf  + '/STEP/Compressor.step')
        print('Pausing 5 seconds for writing Comrpessor STEP.')
        time.sleep(5)

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    def assemble_v8(self, components, *settings):
        # very similar to assemble_v3, except for cylinder intersect
       # Extracting the hub and blade solids from the components
        hub_solid = components[0]  # The hub as a solid
        blade_solids = components[1]  # List of blade solids
        splitter_blade_solids = components[2]  # List of splitter blade solids
        print('hub_solid',hub_solid)
        print('blade_solids',blade_solids)
        print('splitter_blade_solids',splitter_blade_solids)
        # Union the blades with the hub
        for blade in blade_solids:
            hub_solid = hub_solid.union(blade, clean=True, glue=False)
            print('passed a blade union with hub')
        
        # Union the splitter blades with the hub
        for splitter_blade in splitter_blade_solids:
            hub_solid = hub_solid.union(splitter_blade,clean=True, glue=False)
            print('passed a splitter union with hub')

        # Create a larger outer cylinder
        outer_cylinder = cq.Workplane("YZ").circle(2*self.r_4).extrude(self.L_imp)

        # Create a smaller inner cylinder
        inner_cylinder = cq.Workplane("YZ").circle(self.r_4).extrude(self.L_imp)

        # Subtract the smaller cylinder from the larger one to create a hollow cylinder
        hollow_cylinder = outer_cylinder.cut(inner_cylinder, clean=False)

        # Perform the union operation
        hub_solid = hub_solid.cut(hollow_cylinder, clean=False)

        # Union a cylinder (get the orientation right)
        """
        # Create the cylindrical section
        cylinder_height = self.L_imp + self.R_rot/3 #sum([tup[0] for tup in all_coordinates]) #self.L_imp + self.R_rot/3  # Height of the cylinder
        cylinder = (cq.Workplane("YZ")
                    .circle(self.r_4)  # Radius of the cylinder is epsilon
                    .extrude(cylinder_height))
        # Move the cylinder up so that it starts where the hemisphere ends
        cylinder = cylinder.translate((0, 0, 0)).val()

        hub_solid = hub_solid.intersect(cylinder,clean=True)
        """


        """
        hub = (self.sketch_hub
                .revolve(180,(0,0,0),(1,0,0))
                .mirror("XY", union=True)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0), -90)
                .rotate((0,0,0),(0,0,1), 90)) # a Workplane object
        """
        

        # Rest of the assembly code, like saving the combined object
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='Combined Impeller')

        assembly.save(self.cwf  + '/STEP/Compressor.step')
        print('Pausing 5 seconds for writing Comrpessor STEP.')
        time.sleep(5)

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    def assemble_v9(self, components, *settings):
        # very similar to assemble_v3, except for cylinder intersect
       # Extracting the hub and blade solids from the components
        hub_solid = components[0]  # The hub as a solid
        blade_solids = components[1]  # List of blade solids
        splitter_blade_solids = components[2]  # List of splitter blade solids
        print('hub_solid',hub_solid)
        print('blade_solids',blade_solids)
        print('splitter_blade_solids',splitter_blade_solids)
        # Union the blades with the hub
        for blade in blade_solids:
            hub_solid = hub_solid.union(blade, clean=True, glue=False)
            print('passed a blade union with hub')
        
        # Union the splitter blades with the hub
        for splitter_blade in splitter_blade_solids:
            hub_solid = hub_solid.union(splitter_blade,clean=True, glue=False)
            print('passed a splitter union with hub')


        # Create the cylindrical section
        cylinder_height = 2*(self.L_imp + self.R_rot/3) #sum([tup[0] for tup in all_coordinates]) #self.L_imp + self.R_rot/3  # Height of the cylinder
        cylinder = (cq.Workplane("XY")
                    .circle(self.r_4)  # Radius of the cylinder is epsilon
                    .extrude(cylinder_height))
        # Move the cylinder up so that it starts where the hemisphere ends
        cylinder = cylinder.translate((0, 0, -(self.L_imp + self.R_rot/3))).val()

        hub_solid = hub_solid.intersect(cylinder,clean=True, tol=1e-2)


        """
        hub = (self.sketch_hub
                .revolve(180,(0,0,0),(1,0,0))
                .mirror("XY", union=True)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0), -90)
                .rotate((0,0,0),(0,0,1), 90)) # a Workplane object
        """
        

        # Rest of the assembly code, like saving the combined object
        assembly = cq.Assembly(name='Compressor')
        assembly.add(hub_solid, name='Combined Impeller')

        assembly.save(self.cwf  + '/STEP/Compressor.step')
        print('Pausing 5 seconds for writing Comrpessor STEP.')
        time.sleep(5)

        if 'stl' or 'STL' in settings:
            self.helper.convert_step_to_stl(self.cwf + '/STEP/Compressor', self.cwf + '/STL/Compressor')

        return assembly
    
    
    # Defining a method to extract coordinates of the blades from a pickle file
    def blades_coords(self,Element):

        # x,y,z coordinates for the main blade first curve
        self.x_first_curve_bld    = Element['parameters']['comp1']['x_first_curve_bld']
        self.y_first_curve_bld    = Element['parameters']['comp1']['y_first_curve_bld']
        self.z_first_curve_bld    = Element['parameters']['comp1']['z_first_curve_bld']

        # x,y,z coordinates for the main blade first curve
        self.x_second_curve_bld   = Element['parameters']['comp1']['x_second_curve_bld']
        self.y_second_curve_bld   = Element['parameters']['comp1']['y_second_curve_bld']
        self.z_second_curve_bld   = Element['parameters']['comp1']['z_second_curve_bld']

        # x,y,z coordinates for the splitter blade first curve
        self.x_first_curve_split  = Element['parameters']['comp1']['x_first_curve_split']
        self.y_first_curve_split  = Element['parameters']['comp1']['y_first_curve_split']
        self.z_first_curve_split  = Element['parameters']['comp1']['z_first_curve_split']

        # x,y,z coordinates for the splitter blade second curve
        self.x_second_curve_split = Element['parameters']['comp1']['x_second_curve_split']
        self.y_second_curve_split = Element['parameters']['comp1']['y_second_curve_split']
        self.z_second_curve_split = Element['parameters']['comp1']['z_second_curve_split']

        # Storing lists of non-planar blade curves
        
        curve_first_bld     = list(zip(self.x_first_curve_bld,self.y_first_curve_bld,self.z_first_curve_bld))
        curve_second_bld    = list(zip(self.x_second_curve_bld,self.y_second_curve_bld,self.z_second_curve_bld))
        
        curve_first_split   = list(zip(self.x_first_curve_split,self.y_first_curve_split,self.z_first_curve_split))
        curve_second_split  = list(zip(self.x_second_curve_split,self.y_second_curve_split,self.z_second_curve_split))

        # Grouping the curve x,y,z coordinates in tuples   
                
        curve_bld     = [curve_first_bld,curve_second_bld]
        curve_split   = [curve_first_split,curve_second_split]

        return curve_bld, curve_split 