##############################
#importing required libraries#
##############################

import cadquery as cq
import numpy as np
import pandas as pd
from cadquery.selectors import AreaNthSelector
import pickle
from cadquery import *
from cq_warehouse import *
from cq_warehouse.extensions import *

#####################
#loading input files#
#####################

#pickle file with parameters
file = open('Element_23_04_22.pickle','rb')
Element = pickle.load(file)
file.close

#excel file with main blade points
mainblade_file = 'POINT_BLADES1.xls'

#excel file with splitter blade points
splitterblade_file = 'POINT_BLADES2.xls'

#####################################
#split shapes for visualization only#
#####################################

#change boolean to alternate between sides
top_aux=True
bottom_aux=False
top_rotor=True
bottom_rotor=False
top_impeller=True
bottom_impeller=False
top_hub = True
bottom_hub = True

#################################################
#extracting element parameters from pickle in mm#
#################################################

#the coordinate system origin is the center of reference

#diamaters and heights of cylindrical representation
Laenge = [i * 1000 for i in Element['Laenge']]
DI1 = [i * 1000 for i in Element['DI1']]
DI2 = [i * 1000 for i in Element['DI2']]
DI3 = [i * 1000 for i in Element['DI3']]
DA1 = [i * 1000 for i in Element['DA1']]
DA2 = [i * 1000 for i in Element['DA2']]
DA3 = [i * 1000 for i in Element['DA3']]

#tip radius (7mm-35mm)
r_4 = round(Element['parameters']['comp1']['r4'],12)*1000

#tip width (0.1mm-10.5mm)
b_4 = (Element['parameters']['comp1']['b4'])*1000

#inlet hub radius (0.7mm-10.5mm)
r_2h = (Element['parameters']['comp1']['r2h'])*1000

#inlet shoulder radius (0.84mm-24.5mm)
r_2s = (Element['parameters']['comp1']['r2s'])*1000

#inducer inlet radius (0.84mm-35mm)
r_1 = (Element['parameters']['comp1']['r1'])*1000

# diffuser exit radius (7mm-52.5mm)
r_5 = (Element['parameters']['comp1']['r5'])*1000

#blade thickness (0.1mm-0.5mm)
e_bld = (Element['parameters']['comp1']['Blade_e'])*1000

#inducer length (7mm-140mm)
L_ind = (Element['parameters']['comp1']['L_ind'])*1000

#tip clearance (0.001mm-0.158mm)
e_tip = (Element['parameters']['comp1']['Clearance'])*1000

#backface clearance (0.007mm-5.25mm)
e_back = (Element['parameters']['comp1']['Backface'])*1000

#exit blade angle (-45deg-0deg)
beta_4 = Element['parameters']['comp1']['beta4']

if beta_4 ==0:
    beta_4 = 0.01

#number of blades (5blades-11blades)
N_bld = int(Element['parameters']['comp1']['N_blades'])

# radius of rotor
#R_rot = (Element['R_ROT'])*1000

#inlet blade angle (-45deg)
beta_2 = -56

#inlet blade angle shroud (-56deg)
beta_2s = -60

#component types by layer
elem_type1 = Element['elem_type1']
elem_type2 = Element['elem_type2']
elem_type3 = Element['elem_type3']

for k in range(len(Laenge)):
    if elem_type1[k]!='COMP1' or elem_type2[k]!='COMP1' or elem_type3[k]!='COMP1':
            break
# radius of rotor
R_rot = DA3[k]
#R_rot = 5

#defining geometrical parameters
phi = 1.618
b_6 = b_4/phi
a = (r_4/phi)-b_4
b = r_4-r_2s
c = r_4/phi
d = r_4-r_2h
e =(r_4/(phi**2))-b_6
f = r_4-R_rot
L_imp = r_2h+c+b_6+e

#components positions
pos_comp1 = Element['sys_pos']['pos_comp1']


##############################################################
#temporary for manual coupling until we parametrize blades
#manually defining impeller variables (in mm)

# #tip radius (7mm-35mm)
# r_4 = 19

# #inlet shoulder radius (0.84mm-24.5mm)
# r_2s = 2

# #tip width (0.1mm-10.5mm)
# b_4 = 2

# #inducer inlet radius (0.84mm-35mm)
# r_1 = 20

# #inducer hub radius (0.7mm-10.5mm)
# r_2h = 3

# # diffuser exit radius (7mm-52.5mm)
# r_5 = 15

# #blade thickness (0.1mm-0.5mm)
# e_bld = 0.25

# #tip clearance (0.001mm-0.158mm)
# e_tip = 0.01

# #backface clearance (0.007mm-5.25mm)
# e_back = 0.01

# #inducer length (7mm-140mm)
# L_ind = 40

# #exit blade angle (-45deg-0deg)
# beta_4 = -45

# #preventing the twistExtrude() function from failing
# if beta_4 ==0:
#     beta_4 = 0.01

# #inlet blade angle (-45deg)
# beta_2 = -56

# #inlet blade angle shroud (-56deg)
# beta_2s = -60

# #number of blades (5blades-11blades)
# N_bld = 9

# # radius of rotor
# R_rot = 5

# #defining geometrical parameters
# phi = 1.618
# b_6 = b_4/phi
# a = (r_4/phi)-b_4
# b = r_4-r_2s
# c = r_4/phi
# d = r_4-r_2h
# e =(r_4/(phi**2))-b_6
# f = r_4-R_rot
# L_imp = r_2h+c+b_6+e


#defining boolean to either model the rotor or not
show_rotor = True

if show_rotor==True:

    ######################################################
    #automated modeling of cylinder element representation#
    ######################################################

    #initialising disctionaries to store geometrical shapes per function
    layer1={}
    layer2={}
    layer3={}

    #initialising the position at the center of each element along the x axis
    hor_pos=[]
    shift=0
    hor_pos.insert(0,Laenge[0]/2)

    #modeling every element using a loop and the revolve function
    for i in range(len(Laenge)):

        if elem_type1[i]!='COMP1' or elem_type2[i]!='COMP1' or elem_type3[i]!='COMP1':
            #modeling layer 1
            #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                if DI1[i]!=DA1[i]:
                    layer1[i] = (cq.Workplane("XY")
                              .moveTo(hor_pos[i],DI1[i]/2+(DA1[i]-DI1[i])/4)
                              .rect(Laenge[i],(DA1[i]-DI1[i])/2)
                              .revolve(360,(0,0,0),(1,0,0))
                              .split(keepTop=top_aux,keepBottom=bottom_aux)
                              )
            
            #displaying every element
                    # if i in layer1.keys():
                        # show_object(layer1[i],options={"alpha":0, "color": (255, 0, 0)})
        
            #modeling layer 2
            #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                if DI2[i]!=DA2[i]:
        
                    layer2[i] = (cq.Workplane("XY")
                                .moveTo(hor_pos[i],DI2[i]/2+(DA2[i]-DI2[i])/4)
                                .rect(Laenge[i],(DA2[i]-DI2[i])/2)
                                .revolve(360,(0,0,0),(1,0,0))
                                .split(keepTop=top_rotor,keepBottom=bottom_rotor)
                                )
        
            #displaying every element
                    # if i in layer2.keys():
                        # show_object(layer2[i],options={"alpha":0, "color": (0, 255, 0)})
        
            #modeling layer 3
            #modeling every element only if the inner diameter DI is not equal to the outer diameter DA
                if DI3[i]!=DA3[i]:
        
                    layer3[i] = (cq.Workplane("XY")
                                    .moveTo(hor_pos[i],DI3[i]/2+(DA3[i]-DI2[i])/4)
                                    .rect(Laenge[i],(DA3[i]-DI3[i])/2)
                                    .revolve(360,(0,0,0),(1,0,0))
                                    .split(keepTop=top_impeller,keepBottom=bottom_impeller)
                                    )
        
            #displaying every element
                    # if i in layer3.keys():
                        # show_object(layer3[i],options={"alpha":0, "color": (0, 0, 255)})
        else:
            shift += Laenge[i]
    #updating the position at the center of each element along the x axis
        if i<len(Laenge)-1:
            hor_pos.insert(i+1,hor_pos[i]+Laenge[i]/2+Laenge[i+1]/2)





#defining boolean to either model the hub or not
show_hub = True

if show_hub==True:

    #######################################
    #modeling the hub from the pickle file#
    #######################################
    
    #controls the smoothness of the curve
    # step_hub=int(L_imp/b_6)*10
    step_hub=100
    
    #initializing the x and y coordinate arrays
    x_hub = np.linspace(0,L_imp,step_hub)
    y_hub=[]
    y_hub.insert(0,0)
    
    #dividing the hub into 4 parts: hub_1,hub_2,hub_3,hub_4
    #looping over the 4 parts to update the x and y cooridinates
    for i in range(0,len(x_hub)):
    
    #hub_1
        if (x_hub[i]>=0)&(x_hub[i]<=r_2h):
            y_hub[i] = eval('((r_2h**2)-((x_hub[i]-r_2h)**2))**0.5')
            y_hub.insert(i,y_hub[i])
    
    #hub_2
        if (x_hub[i]>=r_2h)&(x_hub[i]<=r_2h+c):
            y_hub[i] = eval('r_4-((d**2)-((d/c)**2)*((x_hub[i]-r_2h)**2))**0.5')
            y_hub.insert(i,y_hub[i])
    
    #hub_3
        if (x_hub[i]>=r_2h+c)&(x_hub[i]<=r_2h+c+b_6):
            y_hub[i]=r_4
            y_hub.insert(i,y_hub[i])
    
    #hub_4
        if (x_hub[i]>=r_2h+c+b_6)&(x_hub[i]<=L_imp):
    
            y_hub[i] = eval('r_4-((f**2)-((f/e)**2)*((x_hub[i]-L_imp)**2))**0.5')
            y_hub.insert(i,y_hub[i])
    
    #grouping the x and y coordinates into tuples
    coords_hub = list(zip(x_hub,y_hub))
    
    # skteching half the hub continuously
    sketch_hub = (cq.Workplane("XY")
                  .polyline(coords_hub,includeCurrent=False)
                  .vLineTo(0)
                  .close()
                  )
    
    #revolving the sketch about the x axis
    hub = (sketch_hub
            .revolve(360,(0,0,0),(1,0,0))
            # .translate((5.25,0,0))
            .split(keepTop=top_hub,keepBottom=bottom_hub)
            )

    #displaying the hub shape
    # show_object(hub,options={"alpha":0, "color": (255,0,0)})

#defining boolean to either model the blades or not
show_blades = False

if show_blades==True:

    ###################################################################
    #retrieving the coordinates of the main blades from the excel file#
    ###################################################################
    
    #storing lists of non-planar blade curves
    points_mainblade = {}
    coords_mainblade={}
    
    #reading the excel file with points of the main blades
    df_mainblade = pd.read_excel(io=mainblade_file,header=None)
    
    #setting index markers to indicate start and end of storing values
    start = None
    end = None
    
    #looping to iterate over all the rows
    for index, row in df_mainblade.iterrows():
    
    #updating index markers every time 'StartCurve' or 'EndCurve' is detected
        if 'StartCurve' in row.values:
            start = index+ 1 #skip the header row
        elif 'EndCurve' in row.values:
            end = index
    
    #storing curve points
            if len(points_mainblade) == 0:
                points_mainblade[0] = df_mainblade.iloc[start:end,:]
            else:
                points_mainblade[len(points_mainblade)] = df_mainblade.iloc[start:end,:]
            start = None
            end = None
    
    #grouping the curve x,y,z coordinates in tuples
    for i in range(len(points_mainblade)):
        coords_mainblade[i] = list(zip(points_mainblade[i].iloc[:,0].astype(float),\
                                        points_mainblade[i].iloc[:,1].astype(float),\
                                        points_mainblade[i].iloc[:,2].astype(float)))
    
    #selecting the coordinates of the last curve
    coords_mainblade_last = list(coords_mainblade.values())[-1]
    
    ##############################################
    #modeling the main blades from the excel file#
    ##############################################
    
    #initializing lists to store sketches and shapes for curves
    sketch_mainblade={}
    
    #Adding the first point of the curve to the end of its point array to close it
    coords_mainblade[0].append(coords_mainblade[0][0])
    coords_mainblade[len(coords_mainblade)-1].append(coords_mainblade[len(coords_mainblade)-1][0])
    
    #initializing lists to store the edges, surfaces, solids, and shells
    mainblade_solid={}
    
    #sketching the first blade curve as an edge
    mainblade_edge0 = cq.Edge.makeSpline([cq.Vector(p) for p in coords_mainblade[0]][0:-1]).close()
    
    #creating a closed surface within the first blade curve edge
    mainblade_face0 = cq.Face.makeNSidedSurface([mainblade_edge0],[])
    
    #sketching the last blade curve as an edge
    mainblade_edge_last = cq.Edge.makeSpline([cq.Vector(p) for p in coords_mainblade[len(coords_mainblade)-1]][0:-1]).close()
    
    #creating a closed surface within the last blade curve edge
    mainblade_face_last = cq.Face.makeNSidedSurface([mainblade_edge_last],[])
    
    #lofting the closed surfaces to produce the overall blade shape
    mainblade_lofted = cq.Solid.makeLoft([mainblade_edge0,mainblade_edge_last],True)
    
    #shelling the faces and edges
    mainblade_shell = cq.Shell.makeShell([mainblade_face0,mainblade_face_last,mainblade_lofted]).fix()
    
    #solidifying the produced shell and rotating
    mainblade_solid[0] = cq.Solid.makeSolid(mainblade_shell).rotate((0,0,0),(0,1,0),90)
    
    #displaying the solidifed main blade
    show_object(mainblade_solid[0])
    
    #looping to model the blade pattern
    for i in range(0,N_bld):
    
    #rotating about the x axis by the corresponding angle
        mainblade_solid[i+1] = (mainblade_solid[i]
                              .transformed((360/N_bld,0,0))
                              )
    
    #displaying the final solid mainblades
        show_object(mainblade_solid[i],options={"alpha":0,"color":(0,0,255)})
    
    #######################################################################
    #retrieving the coordinates of the splitter blades from the excel file#
    #######################################################################
    
    #storing lists of non-planar blade curves
    points_splitterblade = {}
    coords_splitterblade={}
    
    #reading the excel file with points of the splitter blades
    df_splitterblade = pd.read_excel(io=splitterblade_file,header=None)
    
    #setting index markers to indicate start and end of storing values
    start = None
    end = None
    
    #looping to iterate over all the rows
    for index, row in df_splitterblade.iterrows():
    
    #updating index markers every time 'StartCurve' or 'EndCurve' is detected
        if 'StartCurve' in row.values:
            start = index+ 1 #skip the header row
        elif 'EndCurve' in row.values:
            end = index
    
    #storing curve points
            if len(points_splitterblade) == 0:
                points_splitterblade[0] = df_splitterblade.iloc[start:end,:]
            else:
                points_splitterblade[len(points_splitterblade)] = df_splitterblade.iloc[start:end,:]
            start = None
            end = None
    
    #grouping the curve x,y,z coordinates in tuples
    for i in range(len(points_splitterblade)):
        coords_splitterblade[i] = list(zip(points_splitterblade[i].iloc[:,0].astype(float),\
                                        points_splitterblade[i].iloc[:,1].astype(float),\
                                        points_splitterblade[i].iloc[:,2].astype(float)))
    
    #selecting the coordinates of the last curve
    coords_splitterblade_last = list(coords_splitterblade.values())[-1]
    
    ##################################################
    #modeling the splitter blades from the excel file#
    ##################################################
    
    #initializing lists to store sketches and shapes for curves
    sketch_splitterblade={}
    
    #Adding the first point of the curve to the end of its point array to close it
    coords_splitterblade[0].append(coords_splitterblade[0][0])
    coords_splitterblade[len(coords_splitterblade)-1].append(coords_splitterblade[len(coords_splitterblade)-1][0])
    
    #initializing lists to store the edges, surfaces, solids, and shells
    splitterblade_solid={}
    
    #sketching the first blade curve as an edge
    splitterblade_edge0 = cq.Edge.makeSpline([cq.Vector(p) for p in coords_splitterblade[0]][0:-1]).close()
    
    #creating a closed surface within the first blade curve edge
    splitterblade_face0 = cq.Face.makeNSidedSurface([splitterblade_edge0],[])
    
    #sketching the last blade curve as an edge
    splitterblade_edge_last = cq.Edge.makeSpline([cq.Vector(p) for p in coords_splitterblade[len(coords_splitterblade)-1]][0:-1]).close()
    
    #creating a closed surface within the last blade curve edge
    splitterblade_face_last = cq.Face.makeNSidedSurface([splitterblade_edge_last],[])
    
    #lofting the closed surfaces to produce the overall blade shape
    splitterblade_lofted = cq.Solid.makeLoft([splitterblade_edge0,splitterblade_edge_last],True)
    
    #shelling the faces and edges
    splitterblade_shell = cq.Shell.makeShell([splitterblade_face0,splitterblade_face_last,splitterblade_lofted]).fix()
    
    #solidifying the produced shell and rotating
    splitterblade_solid[0] = cq.Solid.makeSolid(splitterblade_shell).rotate((0,0,0),(0,1,0),90)
    
    #displaying the solidifed splitter blade
    show_object(splitterblade_solid[0])
    
    #looping to model the blade pattern
    for i in range(0,N_bld):
    
    #rotating about the x axis by the corresponding angle
        splitterblade_solid[i+1] = (splitterblade_solid[i]
                              .transformed((360/N_bld,0,0))
                              )
    
    #displaying the final solid splitterblades
        show_object(splitterblade_solid[i],options={"alpha":0,"color":(0,0,255)})


#combining the components in the correct position
assembly = cq.Assembly()
assembly.add(hub,color=cq.Color("red"))

for j in range(len(Laenge)):

    if j in layer1.keys():
        assembly.add(layer1[j],loc = cq.Location((abs(L_imp-shift),0,0)),color=cq.Color("red"))

    if j in layer2.keys():
        assembly.add(layer2[j],loc = cq.Location((abs(L_imp-shift),0,0)),color=cq.Color("green"))

    if j in layer3.keys():
        assembly.add(layer3[j],loc = cq.Location((abs(L_imp-shift),0,0)),color=cq.Color("blue"))

# assembly.solve()
show_object(assembly)









