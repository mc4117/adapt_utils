from firedrake import *
from firedrake_adjoint import *
from fenics_adjoint.solving import SolveBlock       # For extracting adjoint solutions
from fenics_adjoint.projection import ProjectBlock  # Exclude projections from tape reading

import datetime
from time import clock

from adapt_utils.adapt.options import DefaultOptions
from adapt_utils.adapt.metric import isotropic_metric, metric_intersection, metric_relaxation


__all__ = ["MeshOptimisation"]


now = datetime.datetime.now()
date = str(now.day) + '-' + str(now.month) + '-' + str(now.year % 2000)


class BaseProblem():
    """
    Base class for solving PDE problems using mesh adaptivity.

    There are three main functionalities:
        * solve PDE;
        * solve adjoint PDE using pyadjoint;
        * adapt mesh based on some error estimator of choice.
    """
    def __init__(self,
                 mesh,
                 finite_element,
                 approach,
                 stab=None,
                 issteady=True,
                 discrete_adjoint=False,
                 op=DefaultOptions(),
                 high_order=False):
        self.mesh = mesh
        self.finite_element = finite_element
        self.approach = approach
        self.stab = stab if stab is not None else 'no'
        self.issteady = issteady
        self.discrete_adjoint = discrete_adjoint
        self.high_order = high_order
        self.op = op
        self.op.approach = approach

        # function spaces and mesh quantities
        self.V = FunctionSpace(self.mesh, self.finite_element)
        self.P0 = FunctionSpace(self.mesh, "DG", 0)
        self.P1 = FunctionSpace(self.mesh, "CG", 1)
        self.P1_vec = VectorFunctionSpace(self.mesh, "CG", 1)
        self.n = FacetNormal(self.mesh)
        self.h = CellSize(self.mesh)

        # prognostic fields
        self.solution = Function(self.V)
        self.adjoint_solution = Function(self.V)

    def set_target_vertices(self, rescaling=0.85, num_vertices=None):
        """
        Set target number of vertices for adapted mesh by scaling the current number of vertices.
        """
        if num_vertices is None:
            num_vertices = self.mesh.num_vertices()
        self.op.target_vertices = num_vertices * rescaling

    def solve(self):
        """
        Solve forward PDE.
        """
        pass


    def objective_functional(self):
        """
        Functional of interest which takes the PDE solution as input.
        """
        pass

    def solve_continuous_adjoint(self):
        """
        Solve the adjoint PDE using a hand-coded continuous adjoint.
        """
        pass

    def solve_discrete_adjoint(self):
        """
        Solve the adjoint PDE in the discrete sense, using pyadjoint.
        """
        if self.issteady:

            # compute some gradient in order to get adjoint solutions
            J = self.objective_functional()
            compute_gradient(J, Control(self.gradient_field))
            tape = get_working_tape()
            solve_blocks = [block for block in tape._blocks if isinstance(block, SolveBlock)
                                                            and not isinstance(block, ProjectBlock)
                                                            and block.adj_sol is not None]
            try:
                assert len(solve_blocks) == 1
            except:
                ValueError("Expected one SolveBlock, but encountered {:d}".format(len(solve_blocks)))

            # extract adjoint solution
            self.adjoint_solution.assign(solve_blocks[0].adj_sol)
            tape.clear_tape()
        else:
            raise NotImplementedError  # TODO

    def solve_adjoint(self):
        """
        Solve adjoint problem using specified method.
        """
        if self.discrete_adjoint:
            self.solve_discrete_adjoint()
        else:
            self.solve_continuous_adjoint()

    def dwp_indication(self):
        """
        Indicate significance by the product of forward and adjoint solutions. This approach was
        used for mesh adaptive tsunami modelling in [Davis and LeVeque, 2016]. Here 'DWP' is used
        to stand for Dual Weighted Primal.
        """
        self.indicator = Function(self.P1)
        self.indicator.project(inner(self.solution, self.adjoint_solution))
        self.indicator.rename('dwp')

    def normalise_indicator(self):
        """
        Given a scalar indicator f and a target number of vertices N, rescale f by
            f := abs(f) * N / norm(f),
        subject to the imposition of minimum and maximum tolerated norms.
        """
        scale_factor = min(max(norm(self.indicator), self.op.min_norm), self.op.max_norm)
        if scale_factor < 1.00001*self.op.min_norm:
            print("WARNING: minimum norm attained")
        elif scale_factor > 0.99999*self.op.max_norm:
            print("WARNING: maximum norm attained")
        self.indicator.interpolate(Constant(self.op.target_vertices/scale_factor)*abs(self.indicator))

    def explicit_estimation(self):
        pass

    def explicit_estimation_adjoint(self):
        pass

    def plot(self):
        """
        Plot current mesh and indicator field, if available.
        """
        di = self.op.directory()
        File(di + 'mesh.pvd').write(self.mesh.coordinates)
        if hasattr(self, 'indicator'):
            name = self.indicator.dat.name
            self.indicator.rename(name + ' indicator')
            File(di + 'indicator.pvd').write(self.indicator)

    def dwr_estimation(self):
        """
        Indicate errors in the objective functional by the Dual Weighted Residual method. This is
        inherently problem-dependent.

        The resulting P0 field should be stored as `self.indicator`.
        """
        pass

    def dwr_estimation_adjoint(self):
        pass

    def get_hessian_metric(self, adjoint=False):
        """
        Compute an appropriate Hessian metric for the problem at hand. This is inherently
        problem-dependent, since the choice of field for adaptation is not universal.

        Hessian metric should be computed and stored as `self.M`.
        """
        pass

    def get_isotropic_metric(self):
        """
        Scale an identity matrix by the indicator field `self.indicator` in order to drive
        isotropic mesh refinement.
        """
        el = self.indicator.ufl_element()
        if (el.family(), el.degree()) != ('Lagrange', 1):
            self.indicator = project(self.indicator, self.P1)
        self.normalise_indicator()
        self.M = isotropic_metric(self.indicator, op=self.op)

    def get_anisotropic_metric(self):
        """
        Apply the approach of [Loseille, Dervieux, Alauzet, 2009] to extract an anisotropic mesh 
        from the Dual Weighted Residual method.
        """
        pass

    def adapt_mesh(self, relaxation_parameter=0.9, prev_metric=None):
        """
        Adapt mesh according to error estimation strategy of choice.
        """
        if self.approach == 'fixed_mesh':
            return
        elif self.approach == 'uniform':
            self.mesh = MeshHierarchy(self.mesh, 1)[0]
            return
        elif self.approach == 'hessian':
            self.get_hessian_metric()
        elif self.approach == 'explicit':
            self.explicit_estimation()
            self.get_isotropic_metric()
        elif self.approach == 'hessian_adjoint':
            self.get_hessian_metric(adjoint=True)
        if self.approach == 'hessian_adjoint':
            self.get_hessian_metric(adjoint=False)
            M = self.M.copy()
            self.get_hessian_metric(adjoint=True)
            self.M = metric_intersection(M, self.M)
        elif self.approach == 'dwp':
            self.dwp_indication()
            self.get_isotropic_metric()
        elif self.approach == 'dwr':
            self.dwr_estimation()
            self.get_isotropic_metric()
        elif self.approach == 'dwr_adjoint':
            self.dwr_estimation_adjoint()
            self.get_isotropic_metric()
        elif self.approach == 'dwr_both':
            self.dwr_estimation()
            self.get_isotropic_metric()
            i = self.indicator.copy()
            self.dwr_estimation_adjoint()
            self.indicator.interpolate(Constant(0.5)*(i+self.indicator))
            self.get_isotropic_metric()
        elif self.approach == 'dwr_averaged':
            self.dwr_estimation()
            self.get_isotropic_metric()
            i = self.indicator.copy()
            self.dwr_estimation_adjoint()
            self.indicator.interpolate(Constant(0.5)*(abs(i)+abs(self.indicator)))
            self.get_isotropic_metric()
        elif self.approach == 'dwr_relaxed':
            self.dwr_estimation()
            self.get_isotropic_metric()
            M = self.M.copy()
            self.dwr_estimation_adjoint()
            self.get_isotropic_metric()
            self.M = metric_relaxation(M, self.M)
        elif self.approach == 'dwr_superposed':
            self.dwr_estimation()
            self.get_isotropic_metric()
            M = self.M.copy()
            self.dwr_estimation_adjoint()
            self.get_isotropic_metric()
            self.M = metric_intersection(M, self.M)
        elif self.approach == 'dwr_anisotropic':
            self.get_anisotropic_metric(adjoint=False)
        elif self.approach == 'dwr_anisotropic_adjoint':
            self.get_anisotropic_metric(adjoint=True)
        elif self.approach == 'dwr_anisotropic_relaxed':
            self.get_anisotropic_metric(adjoint=False)
            M = self.M.copy()
            self.get_anisotropic_metric(adjoint=True)
            self.M = metric_relaxation(M, self.M)
        elif self.approach == 'dwr_anisotropic_superposed':
            self.get_anisotropic_metric(adjoint=False)
            M = self.M.copy()
            self.get_anisotropic_metric(adjoint=True)
            self.M = metric_intersection(M, self.M)
        else:
            raise ValueError("Adaptivity mode {:s} not regcognised.".format(self.approach))

        # Apply metric relaxation, if requested
        self.M_unrelaxed = self.M.copy()
        if prev_metric is not None:
            self.M.project(metric_relaxation(interp(self.mesh, prev_metric), self.M, relaxation_parameter))
        # (Default relaxation of 0.9 following [Power et al 2006])

        # Adapt mesh
        self.mesh = adapt(self.mesh, self.M)

    def interpolate_fields(self):
        """
        Interpolate fields onto the new mesh after a mesh adaptation.
        """
        raise NotImplementedError  # TODO


class MeshOptimisation():
    """
    Loop over all mesh optimisation steps in order to obtain a mesh which is optimal w.r.t. the
    given error estimator for the given PDE problem.
    """
    def __init__(self,
                 problem,
                 mesh,
                 rescaling=0.85,
                 approach='hessian',
                 stab=None,
                 high_order=False,
                 relax=False,
                 outdir='plots/',
                 logmsg='',
                 log=True):

        self.problem = problem
        self.mesh = mesh
        self.rescaling = rescaling
        self.approach = approach
        self.stab = stab if stab is not None else 'no'
        self.high_order = high_order
        self.relax = relax
        self.outdir = outdir
        self.logmsg = logmsg
        self.log = log

        # Default tolerances etc
        self.msg = "Mesh {:2d}: {:7d} cells, objective {:.4e}"
        self.conv_msg = "Converged after {:d} iterations due to {:s}"
        self.maxit = 35
        self.maxit_flag = False
        self.element_rtol = 0.005    # Following [Power et al 2006]
        self.objective_rtol = 0.005

        # Data storage
        self.dat = {'elements': [], 'vertices': [], 'objective': [], 'approach': self.approach}

    def iterate(self):
        M_ = None
        M = None

        # create a log file and spit out parameters used
        if self.log:
            self.logfile = open('{:s}{:s}/optimisation_log'.format(self.outdir, self.approach), 'a+')
            self.logfile.write('\n{:s}{:s}\n\n'.format(date, self.logmsg))
            self.logfile.write('stabilisation: {:s}\n'.format(self.stab))
            self.logfile.write('high_order: {:b}\n'.format(self.high_order))
            self.logfile.write('relax: {:b}\n'.format(self.relax))
            self.logfile.write('maxit: {:d}\n'.format(self.maxit))
            self.logfile.write('element_rtol: {:.3f}\n'.format(self.element_rtol))
            self.logfile.write('objective_rtol: {:.3f}\n\n'.format(self.objective_rtol))

        tstart = clock()
        for i in range(self.maxit):
            print('Solving on mesh {:d}'.format(i))
            tp = self.problem(stab=self.stab,
                              mesh=self.mesh if i == 0 else tp.mesh,
                              approach=self.approach,
                              high_order=self.high_order)

            # Solve
            tp.solve()
            if not self.approach in ('fixed_mesh', 'uniform', 'hessian', 'explicit'):
                tp.solve_adjoint()

            # Extract data
            self.dat['elements'].append(tp.mesh.num_cells())
            self.dat['vertices'].append(tp.mesh.num_vertices())
            self.dat['objective'].append(tp.objective_functional())
            print(self.msg.format(i, self.dat['elements'][i], self.dat['objective'][i]))
            if self.log:
                self.logfile.write('Mesh  {:2d}: elements = {:10d}\n'.format(i, self.dat['elements'][i]))
                self.logfile.write('Mesh  {:2d}: vertices = {:10d}\n'.format(i, self.dat['vertices'][i]))
                self.logfile.write('Mesh  {:2d}:        J = {:.4e}\n'.format(i, self.dat['objective'][i]))

            # Stopping criteria
            if i > 0:
                out = None
                obj_diff = abs(self.dat['objective'][i] - self.dat['objective'][i-1])
                el_diff = abs(self.dat['elements'][i] - self.dat['elements'][i-1])
                if obj_diff < self.objective_rtol*self.dat['objective'][i-1]:
                    out = self.conv_msg.format(i+1, 'convergence in objective functional.')
                elif el_diff < self.element_rtol*self.dat['elements'][i-1]:
                    out = self.conv_msg.format(i+1, 'convergence in mesh element count.')
                elif i >= self.maxit-1:
                    out = self.conv_msg.format(i+1, 'maximum mesh adaptation count reached.')
                    self.maxit_flag = True
                if out is not None:
                    print(out)
                    if self.log:
                        self.logfile.write(out+'\n')
                        tp.plot()
                    break

            # Otherwise, adapt mesh
            tp.set_target_vertices(num_vertices=self.dat['vertices'][0], rescaling=self.rescaling)
            tp.adapt_mesh(prev_metric=M_)
            if self.relax:
                M_ = tp.M_unrelaxed
        self.dat['time'] = clock() - tstart
        print('Time to solution: {:.1f}s'.format(self.dat['time']))
        if self.log:
            self.logfile.close()
