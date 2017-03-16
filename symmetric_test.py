from dolfin import *
from demag_operator import *
import numpy as np

mesh = UnitCubeMesh(5,5,5)

demag_potential = DemagPotentialOperator(mesh)

V = VectorFunctionSpace(mesh, "DG", 0)
v = TestFunction(V)
f = Function(V)
g = Function(V)

A = np.zeros((V.dim(), V.dim()))

for i in range(V.dim()):
  f.vector()[:] = 0.
  f.vector()[i] = 1.
  assemble(inner(f, v)*dx, g.vector())
  u = demag_potential(f)
  A[i,:] = project(- grad(u), V).vector().array()


print A.max()
print (A - A.transpose()).max()

print A[0,:10]
print A[:10,0]
