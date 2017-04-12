from dolfin import *
from dolfin_adjoint import *

import numpy as np
from demag_operator import *

#mesh = Mesh("mesh/cube_pad.xml.gz")
mesh = Mesh("mesh/cubes.xml.gz")

demag = DemagAdjointFunction(mesh, 1, 2, "DG")

m = interpolate(Constant((0,0,1)), demag.function_space)
h_target = Function(demag.function_space, "h_target.xml")
h = demag(m)

dx = Measure('dx', mesh, subdomain_data = MeshFunction('size_t', mesh, 3, mesh.domains()))
J = Functional(0.5 * inner(h - h_target, h - h_target) * dx(2) + 1e-5 * (inner(m, m) - 1.)**2 * dx(1))
m = Control(m)

Jhat = ReducedFunctional(J, m)

m_opt = minimize(Jhat, tol = 1e-12)

File("m.pvd") << m_opt
File("h.pvd") << demag(m_opt)
