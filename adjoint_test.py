from dolfin import *
from dolfin_adjoint import *

import numpy as np
from demag_operator import *
from custom_dolfin_adjoint_function import *

#parameters["adjoint"]["test_derivative"] = True
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


# submeshes and function spaces
mesh_magnet = SubMesh(mesh, 1)
mesh_measure = SubMesh(mesh, 2)
Vmagnet = VectorFunctionSpace(mesh_magnet, "CG", 1)
Vmeasure = VectorFunctionSpace(mesh_measure, "CG", 1)

m = interpolate(Constant((1,0,0)), Vmagnet)

func = CustomDolfinAdjointFunction(forward, adjoint)
h = func(m)

J = Functional(0.5 * inner(h, h) * dx)

# Reduced functional with single control
m = Control(m)

Jhat = ReducedFunctional(J, m)
Jhat.taylor_test(m)# > 1.9
#adj_html("complete.html", "forward")
