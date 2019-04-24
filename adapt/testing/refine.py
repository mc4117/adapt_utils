from firedrake import *
from adapt_utils import *
import numpy as np

# simple mesh of two right angled triangles
mesh = UnitSquareMesh(1, 1)
assert len(mesh.coordinates.dat.data) == 4
P1 = FunctionSpace(mesh, "CG", 1)
P1_ten = TensorFunctionSpace(mesh, "CG", 1)
f = Function(P1)
M = Function(P1_ten)

# check nothing happens when we adapt with the identity metric
M.interpolate(as_matrix([[1/np.sqrt(2), 0], [0, 1/np.sqrt(2)]]))
mesh2 = AnisotropicAdaptation(mesh, M).adapted_mesh
assert np.max(mesh.coordinates.dat.data - mesh2.coordinates.dat.data) < 1e-8

# check isotropic metric does the same thing
f.assign(1/np.sqrt(2))
M = isotropic_metric(f)
mesh2 = AnisotropicAdaptation(mesh, M).adapted_mesh
assert np.max(mesh.coordinates.dat.data - mesh2.coordinates.dat.data) < 1e-8

# check anistropic refinement in x-direction
M2 = anisotropic_refinement(M, direction=0).copy()
mesh2 = AnisotropicAdaptation(mesh, M2).adapted_mesh
assert len(mesh2.coordinates.dat.data) == 6
# TODO: check there are more cells in x-direction

# check anistropic refinement in y-direction
M = isotropic_metric(f)
M3 = anisotropic_refinement(M, direction=1).copy()
mesh2 = AnisotropicAdaptation(mesh, M3).adapted_mesh
assert len(mesh2.coordinates.dat.data) == 6
# TODO: check there are more cells in y-direction

# check metric intersection combines these appropriately
M4 = metric_intersection(M2, M3)
mesh2 = AnisotropicAdaptation(mesh, M4).adapted_mesh
assert len(mesh2.coordinates.dat.data) == 9