from firedrake import *

import numpy as np
import numpy
from numpy import linalg as la
from scipy import linalg as sla

from adapt_utils.options import DefaultOptions
from adapt_utils.adapt.recovery import construct_hessian


__all__ = ["steady_metric", "isotropic_metric", "anisotropic_refinement",
           "metric_intersection", "metric_relaxation", "metric_complexity"]


def steady_metric(f, H=None, mesh=None, noscale=False, op=DefaultOptions()):
    r"""
    Computes the steady metric for mesh adaptation. Based on Nicolas Barral's function
    ``computeSteadyMetric``, from ``adapt.py``, 2016.

    :arg f: P1 solution field.
    :arg H: reconstructed Hessian associated with `f` (if already computed).
    :param op: `Options` class object providing min/max cell size values.
    :return: steady metric associated with Hessian H.
    """
    # NOTE: A P1 field is not actually strictly required
    if H is None:
        H = construct_hessian(f, mesh=mesh, op=op)
    V = H.function_space()
    mesh = V.mesh()
    dim = mesh.topological_dimension()
    assert dim in (2, 3)

    ia2 = 1. / pow(op.max_anisotropy, 2)  # Inverse square max aspect ratio
    ih_min2 = 1. / pow(op.h_min, 2)  # Inverse square minimal side-length
    ih_max2 = 1. / pow(op.h_max, 2)  # Inverse square maximal side-length
    M = Function(V)

    msg = "WARNING: minimum element size reached as {m:.2e}"

    if op.restrict == 'target' or noscale:
        f_min = 1e-6  # Minimum tolerated value for the solution field  # FIXME: seems arbitrary
        if noscale:
            rescale = 1
        else:
            if f is None:
                rescale = op.target
            else:
                rescale = op.target / max(norm(f), f_min)
            #rescale = interpolate(op.target / max_value(abs(f), f_min), P1)

        for i in range(mesh.num_vertices()):

            # Generate local Hessian, avoiding round-off error
            H_loc = rescale*H.dat.data[i]
            mean_diag = 0.5*(H_loc[0][1] + H_loc[1][0])
            H_loc[0][1] = mean_diag
            H_loc[1][0] = mean_diag
            if dim == 3:
                mean_diag = 0.5*(H_loc[0][2] + H_loc[2][0])
                H_loc[0][2] = mean_diag
                H_loc[2][0] = mean_diag
                mean_diag = 0.5*(H_loc[1][2] + H_loc[2][1])
                H_loc[1][2] = mean_diag
                H_loc[2][1] = mean_diag

            # Find eigenpairs and truncate eigenvalues
            lam, v = la.eig(H_loc)
            if dim == 2:
                v1, v2 = v[0], v[1]
            else:
                v1, v2, v3 = v[0], v[1], v[2]
            lam1 = min(ih_min2, max(ih_max2, abs(lam[0])))
            lam2 = min(ih_min2, max(ih_max2, abs(lam[1])))
            if dim == 2:
                lam_max = max(lam1, lam2)
            else:
                lam3 = min(ih_min2, max(ih_max2, abs(lam[2])))
                lam_max = max(lam1, lam2, lam3)
                lam3 = max(lam3, ia2 * lam_max)
            lam1 = max(lam1, ia2 * lam_max)
            lam2 = max(lam2, ia2 * lam_max)
            tol = 0.9999 * ih_min2
            if dim == 2 and (lam[0] >= tol or lam[1] >= tol):
                print(msg.format(m=np.sqrt(min(1./lam[0], 1./lam[1]))))
            elif dim == 3 and (lam[0] >= tol or lam[1] >= tol or lam[2] >= tol):
                print(msg.format(m=np.sqrt(min(1./lam[0], 1./lam[1], 1./lam[2]))))

            # Reconstruct edited Hessian  TODO: condense and generalise notation
            if dim == 2:
                M.dat.data[i][0, 0] = lam1*v1[0]*v1[0] + lam2*v2[0]*v2[0]
                M.dat.data[i][0, 1] = lam1*v1[0]*v1[1] + lam2*v2[0]*v2[1]
                M.dat.data[i][1, 0] = M.dat.data[i][0, 1]
                M.dat.data[i][1, 1] = lam1*v1[1]*v1[1] + lam2*v2[1]*v2[1]
            else:
                M.dat.data[i][0, 0] = lam1*v1[0]*v1[0] + lam2*v2[0]*v2[0] + lam3*v3[0]*v3[0]
                M.dat.data[i][0, 1] = lam1*v1[0]*v1[1] + lam2*v2[0]*v2[1] + lam3*v3[0]*v3[1]
                M.dat.data[i][0, 2] = lam1*v1[0]*v1[2] + lam2*v2[0]*v2[2] + lam3*v3[0]*v3[2]
                M.dat.data[i][1, 0] = M.dat.data[i][0, 1]
                M.dat.data[i][1, 1] = lam1*v1[1]*v1[1] + lam2*v2[1]*v2[1] + lam3*v3[1]*v3[1]
                M.dat.data[i][1, 2] = lam1*v1[1]*v1[2] + lam2*v2[1]*v2[2] + lam3*v3[1]*v3[2]
                M.dat.data[i][2, 0] = M.dat.data[i][0, 2]
                M.dat.data[i][2, 1] = M.dat.data[i][1, 2]
                M.dat.data[i][2, 2] = lam1*v1[2]*v1[2] + lam2*v2[2]*v2[2] + lam3*v3[2]*v3[2]

    elif op.restrict == 'p_norm':
        assert dim == 2  # TODO: 3d
        detH = Function(FunctionSpace(mesh, "CG", 1))

        for i in range(mesh.num_vertices()):

            # Generate local Hessian
            H_loc = H.dat.data[i]
            mean_diag = 0.5 * (H_loc[0][1] + H_loc[1][0])
            H_loc[0][1] = mean_diag
            H_loc[1][0] = mean_diag

            # Find eigenpairs of Hessian and truncate eigenvalues
            lam, v = la.eig(H_loc)
            v1, v2 = v[0], v[1]
            lam1 = max(abs(lam[0]), 1e-10)  # \ To avoid round-off error
            lam2 = max(abs(lam[1]), 1e-10)  # /
            det = lam1 * lam2

            # Reconstruct edited Hessian and rescale
            M.dat.data[i][0, 0] = lam1 * v1[0] * v1[0] + lam2 * v2[0] * v2[0]
            M.dat.data[i][0, 1] = lam1 * v1[0] * v1[1] + lam2 * v2[0] * v2[1]
            M.dat.data[i][1, 0] = M.dat.data[i][0, 1]
            M.dat.data[i][1, 1] = lam1 * v1[1] * v1[1] + lam2 * v2[1] * v2[1]
            M.dat.data[i] *= pow(det, -1. / (2 * op.norm_order + 2))
            detH.dat.data[i] = pow(det, op.norm_order / (2. * op.norm_order + 2))

        # Scale by the target number of vertices and Hessian complexity
        M *= op.target / assemble(detH * dx)

        for i in range(mesh.num_vertices()):
            # Find eigenpairs of metric and truncate eigenvalues
            lam, v = la.eig(M.dat.data[i])
            v1, v2 = v[0], v[1]
            lam1 = min(ih_min2, max(ih_max2, abs(lam[0])))
            lam2 = min(ih_min2, max(ih_max2, abs(lam[1])))
            lam_max = max(lam1, lam2)
            lam1 = max(lam1, ia2 * lam_max)
            lam2 = max(lam2, ia2 * lam_max)
            if (lam[0] >= 0.9999 * ih_min2) or (lam[1] >= 0.9999 * ih_min2):
                print(msg.format(m=np.sqrt(min(1. / lam[0], 1. / lam[1]))))

            # Reconstruct edited Hessian
            M.dat.data[i][0, 0] = lam1 * v1[0] * v1[0] + lam2 * v2[0] * v2[0]
            M.dat.data[i][0, 1] = lam1 * v1[0] * v1[1] + lam2 * v2[0] * v2[1]
            M.dat.data[i][1, 0] = M.dat.data[i][0, 1]
            M.dat.data[i][1, 1] = lam1 * v1[1] * v1[1] + lam2 * v2[1] * v2[1]
    else:
        raise ValueError("Restriction by {:s} not recognised.".format(op.restrict))
    return M


def isotropic_metric(f, bdy=None, noscale=False, op=DefaultOptions()):
    r"""
    Given a scalar error indicator field `f`, construct an associated isotropic metric field.

    :arg f: function to adapt to.
    :param bdy: specify domain boundary to compute metric on.
    :param op: `Options` class providing min/max cell size values.
    :return: isotropic metric corresponding to `f`.
    """
    assert len(f.ufl_element().value_shape()) == 0
    mesh = f.function_space().mesh()
    P1 = FunctionSpace(mesh, "CG", 1)
    dim = mesh.topological_dimension()
    assert dim in (2, 3)

    a2 = pow(op.max_anisotropy, 2)
    h_min2 = pow(op.h_min, 2)
    h_max2 = pow(op.h_max, 2)

    # Scale metric according to restriction strategy
    if noscale or op.restrict == 'p_norm':
        rescale = 1
    else:
        rescale = op.target/min(max(norm(f), op.min_norm), op.max_norm)

    # Project into P1 space and scale
    g = Function(P1)
    g.project(0.5*rescale*f)
    g.interpolate(abs(g))  # ensure non-negative

    # Establish metric
    V = TensorFunctionSpace(mesh, "CG", 1)
    M = Function(V)
    node_set = range(mesh.num_vertices()) if bdy is None else DirichletBC(V, 0, bdy).nodes

    if op.restrict == 'p_norm':
        assert dim == 2  # TODO: 3d case
        detM = Function(P1)
        for i in node_set:
            g.dat.data[i] = max(g.dat.data[i], 1e-8)
            det = g.dat.data[i]**2
            g.dat.data[i] *= pow(det, -1 / 2*op.norm_order + 2)
            detM.dat.data[i] = pow(det, op.norm_order/(2*op.norm_order + 2))
        g *= op.target / assemble(detM*dx)

    for i in node_set:
        alpha = max(1. / h_max2, min(g.dat.data[i], 1. / h_min2))
        alpha = max(alpha, alpha/a2)
        M.dat.data[i][0, 0] = alpha
        M.dat.data[i][1, 1] = alpha
        if dim == 3:
            M.dat.data[i][2, 2] = alpha
        if alpha >= 0.9999 / h_min2:
            print("WARNING: minimum element size reached!")

    return M


# TODO: 3d
def anisotropic_refinement(metric, direction=0):
    r"""
    Anisotropically refine a mesh (or, more precisely, the metric field associated with a mesh)
    in such a way as to approximately half the element size in a canonical direction (x- or y-), by
    scaling of the corresponding eigenvalue.

    :param M: metric to refine.
    :param direction: 0 or 1, corresponding to x- or y-direction, respectively.
    :return: anisotropically refined metric.
    """
    M = metric.copy()
    for i in range(len(M.dat.data)):
        lam, v = la.eig(M.dat.data[i])
        v1, v2 = v[0], v[1]
        lam[direction] *= 4
        M.dat.data[i][0, 0] = lam[0] * v1[0] * v1[0] + lam[1] * v2[0] * v2[0]
        M.dat.data[i][0, 1] = lam[0] * v1[0] * v1[1] + lam[1] * v2[0] * v2[1]
        M.dat.data[i][1, 0] = M.dat.data[i][0, 1]
        M.dat.data[i][1, 1] = lam[0] * v1[1] * v1[1] + lam[1] * v2[1] * v2[1]
    return M


# TODO: 3d
def local_metric_intersection(M1, M2):
    r"""
    Intersect two metrics `M1` and `M2` defined at a particular point in space.
    """
    sqM1 = sla.sqrtm(M1)
    sqiM1 = la.inv(sqM1)  # Note inverse and square root commute whenever both are defined
    lam, v = la.eig(np.dot(np.transpose(sqiM1), np.dot(M2, sqiM1)))
    M12 = np.dot(v, np.dot([[max(lam[0], 1), 0], [0, max(lam[1], 1)]], np.transpose(v)))
    return np.dot(np.transpose(sqM1), np.dot(M12, sqM1))


def metric_intersection(M1, M2, bdy=None):
    r"""
    Intersect a metric field, i.e. intersect (globally) over all local metrics.

    :arg M1: first metric to be intersected.
    :arg M2: second metric to be intersected.
    :param bdy: specify domain boundary to intersect over.
    :return: intersection of metrics M1 and M2.
    """
    V = M1.function_space()
    assert V == M2.function_space()
    M = M1.copy()
    for i in DirichletBC(V, 0, bdy).nodes if bdy is not None else range(V.mesh().num_vertices()):
        M.dat.data[i][:,:] = local_metric_intersection(M1.dat.data[i], M2.dat.data[i])
    return M


def metric_relaxation(M1, M2, alpha=0.5):
    r"""
    Alternatively to intersection, pointwise metric information may be combined using a convex
    combination. Whilst this method does not have as clear an interpretation as metric intersection,
    it has the benefit that the combination may be weighted towards one of the metrics in question.

    :arg M1: first metric to be combined.
    :arg M2: second metric to be combined.
    :param alpha: scalar parameter in [0,1].
    :return: convex combination of metrics M1 and M2 with parameter alpha.
    """
    V = M1.function_space()
    assert V == M2.function_space()
    M = Function(V)
    for i in range(V.mesh().num_vertices()):
        M.dat.data[i][:,:] = alpha*M1.dat.data[i] + (1-alpha)*M2.dat.data[i]
    return M


def metric_complexity(M):
    r"""
    Compute the complexity of a metric, which approximates the number of vertices in a mesh adapted
    based thereupon.
    """
    return assemble(sqrt(det(M))*dx)

