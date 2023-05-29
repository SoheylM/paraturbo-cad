import cadquery as cq
import cq_warehouse.extensions

x_curve = [-1, 0, 1, 1, 1, 0, -1, -1]
y_curve = [2, 1, 1, 1, 1, 1, 2, 2]
z_curve = [3, 3, 2, 1, 0, 1, 2, 3]
coords = list(zip(x_curve, y_curve, z_curve))

edge_1 = cq.Edge.makeSpline([cq.Vector(p) for p in coords][0:-1]).close()

face_1 = cq.Face.makeNSidedSurface([edge_1],[]).thicken(0.5)

CylRadOut1 = 2
CylLen1 = 5
DistCenter1 = 0
L = 1
cylinder1 = cq.Solid.makeCylinder(
     CylRadOut1, 2*CylLen1+2, pnt=cq.Vector(0, DistCenter1+2*L+1, 0), dir=cq.Vector(0, -1, 0)
)

cylinder1=cylinder1.rotate((0,0,0),(1,0,0),90)
cylinder1=cylinder1.transformed((0, 0, 0), (0, 3, 1))
cylinder1=cylinder1.cut(face_1)

show_object(face_1)
show_object(cylinder1)




