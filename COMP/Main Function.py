
import cadquery as cq
from Impeller import *

# from ocp_vscode import show, show_object, reset_show, set_port, set_defaults, get_defaults
# set_port(3939)

# reset_show()
# set_defaults(axes=True, transparent=False, collapse=1, grid=(True, True, True))

Imp = Impeller()

Element = Imp.importpickle('best_solution_Element_2.pickle')
Imp.parameters(Element)
Hub = Imp.hub()

show_object(Hub,options={"alpha":0, "color": (255,0,0)})
