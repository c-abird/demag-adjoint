from dolfin import *
from demag_operator import *

#mesh = Mesh("mesh/cube_pad.xml.gz")
mesh = Mesh("mesh/cubes.xml.gz")
demag_potential = DemagPotentialOperator(mesh)
forward = DirectedDemagFieldOperator(demag_potential, 1, 2)

V = VectorFunctionSpace(mesh, "DG", 0)
m = interpolate(Constant((1, 0, 0)), V)

File("h_target.xml") << forward(m)
File("h_target.pvd") << forward(m)
