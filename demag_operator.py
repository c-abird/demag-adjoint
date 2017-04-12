from dolfin import *
from bempp.api import function_space, GridFunction
from bempp.api.fenics_interface import coupling
from bempp.api.operators.boundary.laplace import double_layer
from bempp.api.operators.boundary.sparse import identity
from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
from bempp.api.common import global_parameters

__all__ = ["DemagPotentialOperator", "DirectedDemagFieldOperator"]

class Struct(object): pass

class DemagPotentialOperator(object):
  def __init__(self, mesh):
    V = FunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    Vv = VectorFunctionSpace(mesh, "DG", 0)
    uv = TrialFunction(Vv)

    self.mesh = mesh
    self._V = V

    # setup u1 problem
    self._u1 = Struct()
    self._u1.result = Function(V)
    A = assemble(inner(grad(u), grad(v)) * dx)
    self._u1.B = assemble(inner(uv, grad(v)) * dx)
    self._u1.solver = LUSolver(A, "mumps")
    self._u1.solver.parameters['reuse_factorization'] = True
    self._u1.solver.parameters['symmetric'] = True

    # setup BEM problem
    self._bem = Struct()
    self._bem.result = Function(V)
    parameters = global_parameters()
    parameters.assembly.boundary_operator_assembly_type = "hmat"
    parameters.quadrature.near.double_order   = 4
    parameters.quadrature.medium.double_order = 3
    parameters.quadrature.far.double_order    = 2
    parameters.hmat.eps = 1e-4

    trace_space, self._bem.trace_matrix = coupling.fenics_to_bempp_trace_data(V)
    dlpOp = double_layer(trace_space, trace_space, trace_space, parameters=parameters, use_projection_spaces=False).weak_form()
    idOp = identity(trace_space, trace_space, trace_space).weak_form()
    invIdOp = InverseSparseDiscreteBoundaryOperator(idOp)
    self._bem.operator = invIdOp * (dlpOp - 0.5 * idOp)

    # setup u2 problem
    self._u2 = Struct()
    self._u2.result = Function(V)
    A = assemble(inner(grad(u), grad(v)) * dx)
    self._u2.b = assemble(Constant(0.) * v * dx)
    # apply fake BCs once to obtain correct system matrix for factorization
    bc = DirichletBC(V, Constant(0.), "on_boundary")
    bc.apply(A)
    self._u2.solver = LUSolver(A, "mumps")
    self._u2.solver.parameters['reuse_factorization'] = True

  def __call__(self, m):
    # u1
    self._u1.solver.solve(self._u1.result.vector(), self._u1.B * m.vector())
    
    # bem
    u2_BEM = self._bem.operator.dot(self._bem.trace_matrix * self._u1.result.vector().array())
    self._bem.result.vector()[:] = self._bem.trace_matrix.transpose() * u2_BEM

    # u2
    bc = DirichletBC(self._V, self._bem.result, "on_boundary")
    bc.apply(self._u2.b)
    self._u2.solver.solve(self._u2.result.vector(), self._u2.b)

    return Function(self._V, self._u1.result.vector() + self._u2.result.vector())

class DirectedDemagFieldOperator(object):
  def __init__(self, demag_operator, source_domain, target_domain):
    self._demag = demag_operator

    mesh = demag_operator.mesh
    V = VectorFunctionSpace(mesh, "DG", 0)
    u = TrialFunction(V)
    v = TestFunction(V)

    cell_domains = MeshFunction('size_t', mesh, 3, mesh.domains())
    dx = Measure('dx', mesh, subdomain_data = cell_domains)

    mass_rec = as_backend_type(assemble(inner(Constant((1,1,1)), v) * dx)).vec()
    mass_rec.reciprocal()
    mass_src = as_backend_type(assemble(inner(Constant((1,1,1)), v) * dx(source_domain))).vec()
    mass_trg = as_backend_type(assemble(inner(Constant((1,1,1)), v) * dx(target_domain))).vec()

    self._Asrc = assemble(inner(u, v) * dx)
    self._Asrc.zero()
    self._Asrc.set_diagonal(PETScVector(mass_rec * mass_src))

    self._Atrg = assemble(inner(u, v) * dx)
    self._Atrg.zero()
    self._Atrg.set_diagonal(PETScVector(mass_rec * mass_trg))
    
    self._V = V

  def __call__(self, m):
    msrc = Function(self._V, self._Asrc * m.vector())
    h = project(- grad(self._demag(msrc)), self._V)
    return Function(self._V, self._Atrg * h.vector())
