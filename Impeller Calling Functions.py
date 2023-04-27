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

Imp.parameters(Element)
Hub = Imp.hub()

Coords_mainblades = Imp.blades_excel('POINT_BLADES1.xls')
Mainblade = Imp.model_blades(Coords_mainblades)
Mainblades = Imp.rotate_blade(Mainblade)

Coords_splitterblades = Imp.blades_excel('POINT_BLADES2.xls')
Splitterblade = Imp.model_blades(Coords_splitterblades)
Splitterblades = Imp.rotate_blade(Splitterblade)

show_object(Hub)
show_object(Mainblades)
show_object(Splitterblades)