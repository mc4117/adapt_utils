from firedrake import *

import numpy
from numpy import linalg as la

from adapt_utils.adapt.metric import steady_metric
from adapt_utils.options import DefaultOptions


__all__ = ["get_eigenpair", "build_metric", "cell_interpolation_error",
           "gradient_interpolation_error", "edge_interpolation_error"]


def get_eigenpair(mesh):
    P0 = FunctionSpace(mesh, "DG", 0)
    P0_vec = VectorFunctionSpace(mesh, "DG", 0)
    P0_ten = TensorFunctionSpace(mesh, "DG", 0)
    lam0 = Function(P0)
    lam1 = Function(P0)
    v0 = Function(P0_vec)
    v1 = Function(P0_vec)
    JJt = Function(P0_ten)
    JJt.interpolate(Jacobian(mesh)*transpose(Jacobian(mesh)))

    for i in range(mesh.num_cells()):
        lam, v = la.eig(JJt.dat.data[i])
        lam0.dat.data[i] = lam[0]
        lam1.dat.data[i] = lam[1]
        v0.dat.data[i][:] = v[0]
        v1.dat.data[i][:] = v[1]

    return [lam0, lam1], [v0, v1]

def build_metric(eigenvalues, eigenvectors, mesh=None):
    if mesh is None:
        mesh = eigenvalues[0].function_space().mesh()
    P0_ten = TensorFunctionSpace(mesh, "DG", 0)
    M = Function(P0_ten)

    for i in range(mesh.num_cells()):
        lam0 = 1/eigenvalues[0].dat.data[i]
        lam1 = 1/eigenvalues[1].dat.data[i]
        v0 = eigenvectors[0].dat.data[i]
        v1 = eigenvectors[1].dat.data[i]
        M.dat.data[i][0, 0] = lam0*v0[0]*v0[0] + lam1*v1[0]*v1[0]
        M.dat.data[i][0, 1] = lam0*v0[0]*v0[1] + lam1*v1[0]*v1[1]
        M.dat.data[i][1, 0] = M.dat.data[i][0, 1]
        M.dat.data[i][1, 1] = lam0*v0[1]*v0[1] + lam1*v1[1]*v1[1]

    return M

def Lij(field, eigenvectors, i, j, mesh=None, op=DefaultOptions()):
    """
    See SIAM paper
    """
    if mesh is None:
        mesh = field.function_space().mesh()
    P0 = FunctionSpace(mesh, "DG", 0)
    P0_ten = TensorFunctionSpace(mesh, "DG", 0)
    I = TestFunction(P0)
    H = project(steady_metric(field, mesh=mesh, noscale=True, op=op), P0_ten)
    triple_product = dot(eigenvectors[i], dot(H, eigenvectors[j]))
    return assemble(I*triple_product*triple_product*dx)

# TODO: test this
def L_matrix(field, eigenvalues, eigenvectors, mesh=None, op=DefaultOptions()):
    if mesh is None:
        mesh = field.function_space().mesh()
    P0_ten = TensorFunctionSpace(mesh, "DG", 0)
    M = Function(P0_ten)
    vals = eigenvalues
    vecs = eigenvectors
    M.interpolate(as_matrix([[vals[0]*vals[0]*Lij(field, vecs, 0, 0, mesh=mesh, op=op),
                               vals[0]*vals[1]*Lij(field, vecs, 0, 1, mesh=mesh, op=op)],
                              [vals[1]*vals[0]*Lij(field, vecs, 1, 0, mesh=mesh, op=op),
                               vals[1]*vals[1]*Lij(field, vecs, 1, 1, mesh=mesh, op=op)]]))
    return M

def cell_interpolation_error(field, eigenvalues, eigenvectors, mesh=None, op=DefaultOptions()):
    if mesh is None:
        mesh = field.function_space().mesh()
    P0 = FunctionSpace(mesh, "DG", 0)
    estimator = Function(P0)
    l = 0
    for i in range(2):
        for j in range(2):
            l = eigenvalues[i]*eigenvalues[j]*Lij(field, eigenvectors, i, j, mesh=mesh, op=op)
    estimator.interpolate(sqrt(l))
    return estimator

def gradient_interpolation_error(field, eigenvalues, eigenvectors, mesh=None, op=DefaultOptions()):
    cell_estimator = cell_interpolation_error(field, eigenvalues, eigenvectors, mesh=mesh, op=op)
    gradient_estimator = Function(cell_estimator.function_space())
    coeff = pow(min_value(eigenvalues[0], eigenvalues[1]), -0.5)
    gradient_estimator.interpolate(coeff*cell_estimator)
    return gradient_estimator

def edge_interpolation_error(field, eigenvalues, eigenvectors, mesh=None, op=DefaultOptions()):
    cell_estimator = cell_interpolation_error(field, eigenvalues, eigenvectors, mesh=mesh, op=op)
    edge_estimator = Function(cell_estimator.function_space())
    coeff = sqrt((eigenvalues[0]+eigenvalues[1])*pow(min_value(eigenvalues[0], eigenvalues[1]), 1.5))
    edge_estimator.interpolate(coeff*cell_estimator)
    return edge_estimator
