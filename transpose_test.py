from dolfin import *
from demag_operator import *
import numpy as np

mesh = Mesh("mesh/cubes.xml.gz")
demag_potential = DemagPotentialOperator(mesh)

# define forward operator
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

submesh1 = SubMesh(mesh, 1)
submesh2 = SubMesh(mesh, 2)
V1 = VectorFunctionSpace(submesh1, "DG", 0)
V2 = VectorFunctionSpace(submesh2, "DG", 0)

#x1 = interpolate(Constant((0,1,0)), V1)
#x2 = interpolate(Constant((0,1,0)), V2)

x1 = Function(V1)
x2 = Function(V2)
x1.vector()[:] = np.random.random(x1.vector().size()) * 1e2
x2.vector()[:] = np.random.random(x2.vector().size()) * 1e2

#print assemble(inner(forward(x1), x2)*dx(submesh2))
#print assemble(inner(backward(x2), x1)*dx(submesh1))

print forward(x1).vector().array().dot(x2.vector().array())
print adjoint(x2).vector().array().dot(x1.vector().array())
