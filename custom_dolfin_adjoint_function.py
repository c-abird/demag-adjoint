from dolfin import *
from dolfin_adjoint import *
import libadjoint

__all__ = ["CustomDolfinAdjointFunction"]

class CustomDolfinAdjointFunction(object):
  def __init__(self, forward, adjoint):
    self._forward = forward
    self._adjoint = adjoint

  def __call__(self, infunc, annotate=None):
    out = self._forward(infunc)

    # Annotate the operation on the dolfin-adjoint tape
    if utils.to_annotate(annotate):
      rhs = CustomDolfinAdjointFunctionRhs(infunc, self._forward, self._adjoint)
      out_dep = adjglobals.adj_variables.next(out)

      solving.register_initial_conditions(zip(rhs.coefficients(), rhs.dependencies()), linear=False)

      if parameters["adjoint"]["record_all"]:
        adjglobals.adjointer.record_variable(out_dep, libadjoint.MemoryStorage(adjlinalg.Vector(out)))

      identity = utils.get_identity_block(out.function_space())
      eq = libadjoint.Equation(out_dep, blocks=[identity], targets=[out_dep], rhs=rhs)
      cs = adjglobals.adjointer.register_equation(eq)

      solving.do_checkpoint(cs, out_dep, rhs)

    return out

class CustomDolfinAdjointFunctionRhs(libadjoint.RHS):
  def __init__(self, infunc, forward, adjoint):
    self.infunc = infunc
    self._forward = forward
    self._adjoint = adjoint
    self.deps = [adjglobals.adj_variables[self.infunc]]

  def __call__(self, dependencies, values):
    out = self._forward(values[0].data)
    return adjlinalg.Vector(out)

  def derivative_action(self, dependencies, values, variable, contraction_vector, hermitian):
    if hermitian:
      out = self._adjoint(contraction_vector.data)
    else:
      raise Exception()
      out = self._adjoint(contraction_vector.data)

    return adjlinalg.Vector(out)

  def dependencies(self):
    # What does this custom function depend on?
    return self.deps 

  def coefficients(self):
    return [self.infunc]

  def __str__(self):
    return "CustomDolfinAdjointFunction"
