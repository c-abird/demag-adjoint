from dolfin import *
from bempp.api import function_space, GridFunction
from bempp.api.fenics_interface import coupling
from bempp.api.operators.boundary.laplace import double_layer
from bempp.api.operators.boundary.sparse import identity
from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
from bempp.api.common import global_parameters
from custom_dolfin_adjoint_function import *

__all__ = ["DemagPotentialOperator", "DirectedDemagFieldOperator", "DemagAdjointFunction"]

class Struct(object): pass

class DemagPotentialOperator(object):
  def __init__(self, mesh, function_space = "DG"):
    V = FunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    if function_space == "DG":
      Vv = VectorFunctionSpace(mesh, "DG", 0)
    elif function_space == "CG":
      Vv = VectorFunctionSpace(mesh, "CG", 1)
    else:
      raise Exception()
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
  def __init__(self, demag_operator, source_domain, target_domain, function_space = "DG"):
    self._demag = demag_operator

    mesh = demag_operator.mesh
    if function_space == "DG":
      V = VectorFunctionSpace(mesh, "DG", 0)
    elif function_space == "CG":
      V = VectorFunctionSpace(mesh, "CG", 1)
    else:
      raise Exception()

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

class DemagAdjointFunction(CustomDolfinAdjointFunction):
  def __init__(self, mesh, source_domain, target_domain, function_space = "DG"):
    if function_space == "DG":
      self.function_space = VectorFunctionSpace(mesh, "DG", 0)
    elif function_space == "CG":
      self.function_space = VectorFunctionSpace(mesh, "CG", 1)
    else:
      raise Exception()

    demag_potential = DemagPotentialOperator(mesh, function_space)
    forward = DirectedDemagFieldOperator(demag_potential, source_domain, target_domain, function_space)

    # define adjoint operator
    backward = DirectedDemagFieldOperator(demag_potential, target_domain, source_domain, function_space)
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


    super(DemagAdjointFunction, self).__init__(forward, adjoint)
