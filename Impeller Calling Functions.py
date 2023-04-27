import cadquery as cq
import os
cwf = os.getcwd()

from ocp_vscode import show, show_object, reset_show, set_port, set_defaults, get_defaults
set_port(3939)

reset_show() # use for reapeated shift-enter execution to clean object buffer
set_defaults(axes=True, transparent=False, collapse=1, grid=(True, True, True))

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

#Help.assemble((Hub,Mainblades,Splitterblades))

# show_object(Hub)
# show_object(Mainblades)
# show_object(Splitterblades)
