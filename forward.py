from dolfin import *
from demag_operator import *

#mesh = Mesh("mesh/cube_pad.xml.gz")
mesh = Mesh("mesh/cubes.xml.gz")
demag = DemagAdjointFunction(mesh, 1, 2, "DG")
m = interpolate(Constant((1, 0, 0)), demag.function_space)

File("h_target.xml") << demag(m)
File("h_target.pvd") << demag(m)
