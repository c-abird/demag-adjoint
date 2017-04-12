from dolfin import *
from dolfin_adjoint import *

import numpy as np
from demag_operator import *
from custom_dolfin_adjoint_function import *

#mesh = Mesh("mesh/cube_pad.xml.gz")
mesh = Mesh("mesh/cubes.xml.gz")

demag_potential = DemagPotentialOperator(mesh)
forward = DirectedDemagFieldOperator(demag_potential, 1, 2)

# define adjoint operator
backward = DirectedDemagFieldOperator(demag_potential, 2, 1)
def adjoint(m):
  u = TrialFunction(m.function_space())
  v = TestFunction(m.function_space())
  inp = Function(m.function_space())
  A = assemble(inner(u, v)*dx)
  solve(A, inp.vector(), m.vector())

  h = backward(inp)

  u = TrialFunction(h.function_space())
  v = TestFunction(h.function_space())
  result = Function(h.function_space())
  assemble(inner(h, v)*dx, result.vector())

  return result

V = VectorFunctionSpace(mesh, "DG", 0)
dx = Measure('dx', mesh, subdomain_data = MeshFunction('size_t', mesh, 3, mesh.domains()))

m = interpolate(Constant((0,0,1)), V)

demag = CustomDolfinAdjointFunction(forward, adjoint)

h_target = Function(V, "h_target.xml")
#h_target = Constant((0,0,1))
h = demag(m)

J = Functional(0.5 * inner(h - h_target, h - h_target) * dx(2) + 1e-20 * inner(m, m) * dx(1))
m = Control(m)

Jhat = ReducedFunctional(J, m)

Jhat.taylor_test(m)
exit()
#m_opt = minimize(Jhat, tol = 1e-20)

File("m.pvd") << m_opt
File("h.pvd") << demag(m_opt)
