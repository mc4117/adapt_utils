from firedrake import *

from adapt_utils.adapt.kernels import *


__all__ = ["AdaptiveMesh"]


class AdaptiveMesh():
    """
    Wrapper which adds extra features to mesh.
    """
    def __init__(self, mesh, levels=0):
        """
        `AdaptMesh` object is initialised as the basis of a `MeshHierarchy`.
        """
        self.levels = levels
        self.hierarchy = MeshHierarchy(mesh, levels)
        self.mesh = self.hierarchy[0]
        self.dim = self.mesh.topological_dimension()
        assert self.dim in (2, 3)
        if levels > 0:
            self.refined_mesh = self.hierarchy[1]

        self.normal = FacetNormal(self.mesh)
        if self.dim == 2:
            self.tangent = as_vector([-self.normal[1], self.normal[0]])  # Tangent vector
        elif self.dim == 3:
            raise NotImplementedError  # TODO: Get a tangent vector in 3D
        else:
            raise NotImplementedError
        self.facet_area = FacetArea(self.mesh)
        self.h = CellSize(self.mesh)

    def save_plex(self, filename):
        """
        Save mesh in DMPlex format.
        """
        viewer = PETSc.Viewer().createHDF5(filename, 'r')
        viewer(self.mesh._plex)

    def load_plex(self, filename):
        """
        Load mesh from DMPlex format. The `MeshHierarchy` is reinstated.
        """
        newplex = PETSc.DMPlex().create()
        newplex.createFromFile(filename)
        self.__init__(Mesh(newplex), levels=self.levels)

    def adapt(self, metric):
        """
        Adapt mesh using a specified metric. The `MeshHierarchy` is reinstated.
        """
        self.__init__(adapt(self.mesh, metric), levels=self.levels)

    def copy(self):
        return AdaptiveMesh(Mesh(Function(self.mesh.coordinates)), levels=self.levels)

    def get_edge_lengths(self):
        """
        For each element, find the lengths of associated edges, stored in a HDiv trace field.

        NOTE: The plus sign is arbitrary and could equally well be chosen as minus.
        """
        HDivTrace = FunctionSpace(self.mesh, "HDiv Trace", 0)
        v, u = TestFunction(HDivTrace), TrialFunction(HDivTrace)
        self.edge_lengths = Function(HDivTrace, name="Edge lengths")
        mass_term = v('+')*u('+')*dS + v*u*ds
        rhs = v('+')*self.facet_area*dS + v*self.facet_area*ds
        solve(mass_term == rhs, self.edge_lengths)

    def get_edge_vectors(self):
        """
        For each element, find associated edge vectors, stored in a HDiv trace field.

        NOTES:
          * The plus sign is arbitrary and could equally well be chosen as minus.
          * The sign of the returned vectors is arbitrary and could equally well take the minus sign.
        """
        HDivTrace_vec = VectorFunctionSpace(mesh, "HDiv Trace", 0)
        v, u = TestFunction(HDivTrace_vec), TrialFunction(HDivTrace_vec)
        self.edge_vectors = Function(HDivTrace_vec, name="Edge vectors")
        mass_term = inner(v('+'), u('+'))*dS + inner(v, u)*ds
        rhs = inner(v('+'), self.tangent('+')*self.facet_area)*dS
        rhs += inner(v, self.tangent*self.facet_area)*ds
        solve(mass_term == rhs, self.edge_vectors)

    def get_maximum_length_edge(self):
        """
        For each element, find the associated edge of maximum length.
        """
        self.get_edge_lengths()
        self.get_edge_vectors()
        P0_vec = VectorFunctionSpace(self.mesh, "DG", 0)
        self.maximum_length_edge = Function(P0_vec, name="Maximum length edge")
        par_loop(get_maximum_length_edge(self.dim), dx, {'edges': (self.edge_lengths, READ),
                                                         'vectors': (self.edge_vectors, READ),
                                                         'max_vector': (self.maximum_length_edge, RW)
                                                        })

    def get_cell_metric(self):
        """
        Compute cell metric associated with mesh.

        Based on code by Lawrence Mitchell.
        """
        P0_ten = TensorFunctionSpace(self.mesh, "DG", 0)
        J = interpolate(Jacobian(self.mesh), P0_ten)
        self.cell_metric = Function(P0_ten, name="Cell metric")
        kernel = eigen_kernel(singular_value_decomposition, dim)
        op2.par_loop(kernel, P0_ten.node_set, self.cell_metric.dat(op2.INC), J.dat(op2.READ))

    def anisotropic_h(self, u):
        """
        Measure of element size recommended in [Nguyen et al., 2009]: maximum edge length, projected
        onto the velocity field `u`.
        """
        try:
            assert isinstance(u, Constant) or isinstance(u, Function)
        except AssertionError:
            raise ValueError("Velocity field should be either `Function` or `Constant`.")
        self.get_maximum_length_edge()
        v = self.maximum_length_edge
        P0 = FunctionSpace(self.mesh, "DG", 0)
        h = interpolate((u[0]*v[0] + u[1]*v[1])/sqrt(dot(u, u)), P0)
        return h.vector().gather().max()  # TODO: Spatially varying version

