import cadquery as cq
import os
cwf = os.getcwd()

import sys
sys.path.append(cwf  + '/COMP')

from Impeller import *
from Helper import *

Help = Helper()

Element = Help.importpickle('Element_23_04_22.pickle')

Imp = Impeller()

Imp.parameters_rotor(Element)
#Imp.parameters_impeller(Element)


#tip radius (7mm-35mm)
r_4 = 19
#inlet shoulder radius (0.84mm-24.5mm)
r_2s = 2
#tip width (0.1mm-10.5mm)
b_4 = 2
#inducer inlet radius (0.84mm-35mm)
r_1 = 20
#inducer hub radius (0.7mm-10.5mm)
r_2h = 3
# diffuser exit radius (7mm-52.5mm)
r_5 = 15
#blade thickness (0.1mm-0.5mm)
e_bld = 0.25
#tip clearance (0.001mm-0.158mm)
e_tip = 0.01
#backface clearance (0.007mm-5.25mm)
e_back = 0.01
#inducer length (7mm-140mm)
L_ind = 40
#exit blade angle (-45deg-0deg)
beta_4 = -45
#inlet blade angle (constant)
beta_2 = -56
#inlet blade angle shroud (constant)
beta_2s = -60
#number of blades (5blades-11blades)
N_bld = 9
# radius of rotor (mm)
R_rot = 5


Imp.manualparams_impeller(r_4,r_2s,beta_4,b_4,r_1,r_2h,r_5,e_bld,e_tip,e_back,L_ind,beta_2,beta_2s,N_bld,R_rot)
Hub = Imp.hub()
Imp.settings_hub(True,True)

Coords_mainblades = Imp.blades_excel('POINT_BLADES1.xls')
Mainblade = Imp.model_blades(Coords_mainblades)
Mainblades = Imp.rotate_blade(Mainblade)

Coords_splitterblades = Imp.blades_excel('POINT_BLADES2.xls')
Splitterblade = Imp.model_blades(Coords_splitterblades)
Splitterblades = Imp.rotate_blade(Splitterblade)

#Help.assemble((Hub,Mainblades,Splitterblades))

show_object(Hub)
show_object(Mainblades)
show_object(Splitterblades)
