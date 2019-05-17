from firedrake import *

from adapt_utils.options import DefaultOptions


__all__ = ["iso_P2", "multi_adapt"]


def iso_P2(mesh):
    r"""
    Uniformly refine a mesh (in each canonical direction) using an iso-P2 refinement. That is, nodes
    of a quadratic element on the initial mesh become vertices of the new mesh.
    """
    return MeshHierarchy(mesh, 1).__getitem__(1)


def multi_adapt(metric, op=DefaultOptions()):
    r"""
    Adapt mesh multiple times, by repeatedly projecting the metric into the new space.

    This should be done more than once, but at most four times. The more steps applied, the larger
    the errors attributed to the projection.
    """
    for i in range(op.num_adapt):
        mesh = metric.function_space().mesh()
        newmesh = adapt(mesh, metric)
        if i < op.num_adapt-1:
            V = TensorFunctionSpace(newmesh, "CG", 1)
            metric = project(metric, V)
    return newmesh
