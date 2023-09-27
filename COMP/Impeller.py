# Importing required libraries
import cadquery as cq
import numpy as np
import pandas as pd
#from cadquery import *
#from cq_warehouse import *
#from cq_warehouse.extensions import *
import os

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

        # NO NEEDED
        """
        r_2h = self.r_2h
        r_4 = self.r_4
        b_4 = self.b_4
        R_rot = self.R_rot  
        
        phi = 1.618
        self.b_6 = b_4/phi
        self.c = r_4/phi
        self.e =(r_4/(phi**2))-self.b_6
        self.L_imp = r_2h+self.c+self.b_6+self.e
        """

        # Tolerance to control the smoothness of the hub profile
        self.delta_x = self.b_6/3
        #print('b_6',self.b_6)

        # BUGGY BECAUSE delta_x NOT CONSTANT BETWEEN x_1[-1] and x_2[0], x_2[-1] and x_3[0] etc.
        """
        # Define the intervals
        x_1 = np.arange(0, self.r_2h, self.delta_x)
        print('x_1',x_1)
        x_2 = np.arange(self.r_2h, self.r_2h + self.c, self.delta_x)
        print('x_2',x_2)
        x_3 = np.arange(self.r_2h + self.c, self.r_2h + self.c + self.b_6, self.delta_x)
        print('x_3',x_3)
        x_4 = np.arange(self.r_2h + self.c + self.b_6, self.L_imp, self.delta_x)
        print('x_4',x_4)
        x_5 = np.arange(self.L_imp, self.L_imp + self.R_rot/3, self.delta_x)
        print('x_5',x_5)
        """

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
        
        # working but slight bug with disconnected blade which I suspect to be linked to the addition of zero at beginnig of y_hub
        """

        # Initializing the x and y coordinate arrays
        n_hub = int(np.round((L_imp+(R_rot/3))/delta_x))
        self.x_hub = np.linspace(0,L_imp+(R_rot/3), num = n_hub)
        self.y_hub=[]
        self.y_hub.insert(0,0)

        # Initialize the three lists to mix spline and polyline
        self.coords_hub_before_3 = []
        self.coords_hub_3 = []
        self.coords_hub_after_3 = []

        # Dividing the hub into 5 parts: hub_1,hub_2,hub_3,hub_4,hub_5
        # Looping over the 5 parts to update the x and y coordinates
        for i in range(0,len(self.x_hub)):

        # Hub_1
            if (self.x_hub[i]>=0)&(self.x_hub[i]<r_2h):
                self.y_hub[i] = eval('((self.r_2h**2)-((self.x_hub[i]-self.r_2h)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])
                self.coords_hub_before_3.append((self.x_hub[i], self.y_hub[i]))
                print('self.x_hub[i] hub1',self.x_hub[i])
                print('self.y_hub[i] hub1',self.y_hub[i])

        # Hub_2
            if (self.x_hub[i]>=r_2h)&(self.x_hub[i]<r_2h+c):
                self.y_hub[i] = eval('self.r_4-((self.d**2)-((self.d/self.c)**2)*((self.x_hub[i]-self.r_2h)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])
                self.coords_hub_before_3.append((self.x_hub[i], self.y_hub[i]))
                print('self.x_hub[i] hub2',self.x_hub[i])
                print('self.y_hub[i] hub2',self.y_hub[i])

        # Hub_3
            if (self.x_hub[i]>=r_2h+c)&(self.x_hub[i]<r_2h+c+b_6):
                self.y_hub[i] = r_4
                self.y_hub.insert(i,self.y_hub[i])
                self.coords_hub_3.append((self.x_hub[i], self.y_hub[i]))
                print('self.x_hub[i] hub3',self.x_hub[i])
                print('self.y_hub[i] hub3',self.y_hub[i])
                

        # Hub_4
            if (self.x_hub[i]>=r_2h+c+b_6)&(self.x_hub[i]<L_imp):
                self.y_hub[i] = eval('self.r_4-((self.f**2)-((self.f/self.e)**2)*((self.x_hub[i]-self.L_imp)**2))**0.5')
                self.y_hub.insert(i,self.y_hub[i])
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))
                print('self.x_hub[i] hub4',self.x_hub[i])
                print('self.y_hub[i] hub4',self.y_hub[i])

        # Hub_5
            if (self.x_hub[i]>=L_imp)&(self.x_hub[i]<=L_imp+(R_rot/3)):
                self.y_hub[i] = R_rot
                self.y_hub.insert(i,self.y_hub[i])
                self.coords_hub_after_3.append((self.x_hub[i], self.y_hub[i]))
                print('self.x_hub[i] hub5',self.x_hub[i])
                print('self.y_hub[i] hub5',self.y_hub[i])


        # Grouping the x and y coordinates into tuples
        self.coords_hub = list(zip(self.x_hub,self.y_hub))
        """

        # Create the profile with spline and polyline
        self.sketch_hub = (cq.Workplane("XY")
            .spline(self.coords_hub_before_3, includeCurrent=False)
            .polyline(self.coords_hub_3, includeCurrent=True)
            .spline(self.coords_hub_after_3, includeCurrent=True)
            .vLineTo(0)
            .close())

        """
        # Skteching half the hub continuously
        self.sketch_hub = (cq.Workplane("XY")
                    .polyline(self.coords_hub,includeCurrent=False)
                    .vLineTo(0)
                    .close())
        """

        # Calculating the shift in the position of the true impeller model
        self.shift = 0
        for i in range(len(self.Laenge)):
            if self.elem_type1[i]=='COMP1' or self.elem_type2[i]=='COMP1' or self.elem_type3[i]=='COMP1':
                self.shift += self.Laenge[i]

        # Revolving the sketch about the rotation axis
        hub = (self.sketch_hub
                .revolve(360,(0,0,0),(1,0,0))
                .split(keepTop=True,keepBottom=bottom_hub)
                .translate((-abs(self.L_imp-self.shift),0,0))
                .rotate((0,0,0),(0,1,0),-90)
                .rotate((0,0,0),(0,0,1),90))
        
        #print('self.L_imp-self.shift',self.L_imp-self.shift)

        assembly = cq.Assembly(name = 'Turboompressor Hub')
        assembly.add(hub,color=cq.Color('red3'), name = 'Hub')

        return assembly
    
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
        blade_solid[0] = cq.Solid.makeSolid(blade_shell).translate((0,0,-abs(self.L_imp-self.shift)))
        #print('abs(self.L_imp-self.shift)',abs(self.L_imp-self.shift))

        return blade_solid
    
    def trim_blades(self, blade):

        # Create a larger outer cylinder
        outer_cylinder = cq.Solid.makeCylinder(2*self.r_4, self.L_imp, cq.Vector(0, 0, -abs(self.L_imp-self.shift)))

        # Create a smaller inner cylinder
        inner_cylinder = cq.Solid.makeCylinder(self.r_4, self.L_imp, cq.Vector(0, 0, -abs(self.L_imp-self.shift)))

        # Subtract the smaller cylinder from the larger one to create a hollow cylinder
        hollow_cylinder = outer_cylinder.cut(inner_cylinder)

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