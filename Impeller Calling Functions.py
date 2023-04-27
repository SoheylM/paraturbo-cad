
import cadquery as cq
from Impeller import *

from ocp_vscode import show, show_object, reset_show, set_port, set_defaults, get_defaults
set_port(3939)

reset_show()
set_defaults(axes=True, transparent=False, collapse=1, grid=(True, True, True))

Imp = Impeller()

Element = Imp.importpickle('Element_23_04_22.pickle')
Imp.parameters(Element)
Hub = Imp.hub()

show_object(Hub)



