from firedrake import *
from firedrake_adjoint import *
from firedrake.slate.slac.compiler import PETSC_DIR
import firedrake.dmplex as dmplex
import numpy as np


__all__ = ["anisotropic_stabilisation"]


def cell_metric(mesh, metric=None):
    """
    Compute cell metric associated with mesh.

    Based on code by Lawrence Mitchell.
    """
    #print("Making cell metric on %s" % mesh)
    dim = mesh.topological_dimension()
    assert dim in (2, 3)
    P0_ten = TensorFunctionSpace(mesh, "DG", 0)
    J = interpolate(Jacobian(mesh), P0_ten)
    metric = metric or Function(P0_ten, name="CellMetric")
    d = str(dim)
    kernel_str = """
#include <Eigen/Dense>

void polar(double A_[%s], const double * B_) {
  Eigen::Map<Eigen::Matrix<double, %s, %s, Eigen::RowMajor> > A((double *)A_);
  Eigen::Map<Eigen::Matrix<double, %s, %s, Eigen::RowMajor> > B((double *)B_);
  Eigen::JacobiSVD<Eigen::Matrix<double, %s, %s, Eigen::RowMajor> > svd(B, Eigen::ComputeFullV);

  A += svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
}""" % (str(dim*dim), d, d, d, d, d, d)
    kernel = op2.Kernel(kernel_str, "polar", cpp=True, include_dirs=["%s/include/eigen3" % d for d in PETSC_DIR])
    op2.par_loop(kernel, P0_ten.node_set, metric.dat(op2.INC), J.dat(op2.READ))
    return metric


def anisotropic_stabilisation(u, mesh=None):
    """
    Compute anisotropic stabilisation coefficient using `cell_metric` and the velocity field.
    """
    if mesh is None:
        mesh = u.function_space().mesh()
    M = cell_metric(mesh)
    P0 = FunctionSpace(mesh, "DG", 0)
    h = dot(u, dot(M, u))  # chosen measure of cell size which accounts for anisotropy
    tau = Function(P0)
    tau.interpolate(0.5*h/sqrt(inner(u, u)))
    return tau


def supg_coefficient(u, nu, mesh=None, anisotropic=False):
    r"""
    Compute SUPG stabilisation coefficent for an advection diffusion problem. There are two modes
    in which this can be calculated, as determined by the Boolean parameter `anisotropic`:

    In isotropic mode, we use the cell diameter as our measure of element size :math:`h_K`.

    In anisotropic mode, we loop over each element of the mesh, project the edge of maximal length
    into a vector space spanning the velocity field `u` and take the length of this projected
    edge as the measure of element size.

    In both cases, we compute the stabilisation coefficent as

..  math::
    \tau = \frac{h_K}{2\|\textbf{u}\|}

    :arg u: velocity field associated with advection equation being solved.
    :arg nu: diffusivity of fluid.
    :kwarg mesh: mesh upon which problem is defined.
    :kwarg anisotropic: toggle between isotropic and anisotropic mode.
    """
    if mesh is None:
        mesh = u.function_space().mesh()
    h = anisotropic_h(u, mesh) if anisotropic else CellSize(mesh)
    Pe = 0.5*sqrt(inner(u, u))*h/nu
    tau = 0.5*h/sqrt(inner(u, u))
    return tau*min_value(1, Pe/3)


def anisotropic_h(u, mesh=None):
    if mesh is None:
        mesh = u.function_space().mesh()
    plex = mesh._plex
    func = isinstance(u, Function)  # determine if u is a Function or a Constant

    # create section describing global numbering of vertices
    dim = mesh.topological_dimension()
    entity_dofs = np.zeros(dim+1, dtype=np.int32)
    entity_dofs[0] = mesh.geometric_dimension()
    section = dmplex.create_section(mesh, entity_dofs)

    # extract strata from plex
    v_start, v_end = plex.getHeightStratum(2)  # vertices
    f_start, f_end = plex.getHeightStratum(1)  # edges/facets
    c_start, c_end = plex.getHeightStratum(0)  # cells/elements

    dat = {}
    for c in range(c_start, c_end):                           # loop over elements
        dat[c] = {}
        max_norm = 0.
        edges = plex.getCone(c)                               # get edges of element
        for f in edges:                                       # loop over edge of element
            assert f >= f_start and f < f_end
            vertices = plex.getCone(f)                        # get vertices of edge
            crds = []
            if func:
                func_vals = []
            for v in vertices:                                # loop over vertices of edge
                assert v >= v_start and v < v_end
                off = section.getOffset(v)//2
                crds.append(mesh.coordinates.dat.data[off])   # get coordinates of vertex
                if func:
                    func_vals.append(u.dat.data[off])         # get corresp. values of u
            vec = crds[1] - crds[0]                           # get vector describing edge
            nrm = np.sqrt(np.dot(vec, vec))                   # take norm of vector
            if func:
                avg_func = (0.5*func_vals[0] + func_vals[1])  # take average along edge

            if nrm >= max_norm:                               # find edge with max length
                dat[c]['edge'] = vec
                if func:
                    func_val = avg_func
        dat[c]['u'] = func_val if func else u.dat.data

    # create section describing global numbering of elements
    entity_dofs = np.zeros(dim+1, dtype=np.int32)
    entity_dofs[2] = mesh.geometric_dimension()
    section = dmplex.create_section(mesh, entity_dofs)

    # project edge with max length into vector space spanned by velocity field
    P0 = FunctionSpace(mesh, "DG", 0)
    h = Function(P0)
    for c in dat.keys():
        e = dat[c]['edge']
        u_e = dat[c]['u']
        projected_edge = np.array([e[0]*abs(u_e[0]), e[1]*abs(u_e[1])])
        idx = section.getOffset(c)//2
        h.dat.data[idx] = np.sqrt(np.dot(projected_edge, projected_edge))

    return h