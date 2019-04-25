from firedrake import *
from firedrake_adjoint import *
from thetis.configuration import *
import math
from scipy.special import kn

from adapt_utils.options import Options
from adapt_utils.misc import *


__all__ = ["TracerOptions", "PowerOptions", "TelemacOptions", "TelemacOptions_Centred",
           "LeVequeOptions"]


class TracerOptions(Options):
    """
    Default parameter class for TracerProblems.
    """

    # Domain
    nx = PositiveInteger(4, help="Mesh resolution in x- and y-directions.").tag(config=True)
    target_vertices = PositiveFloat(1000., help="Target number of vertices (not an integer!)")

    # Timestepping
    dt = PositiveFloat(0.1, help="Timestep").tag(config=True)
    start_time = NonNegativeFloat(0., help="Start of time window of interest").tag(config=True)
    end_time = PositiveFloat(60., help="End of time window of interest").tag(config=True)
    dt_per_export = PositiveFloat(10, help="Number of timesteps per export").tag(config=True)
    dt_per_remesh = PositiveFloat(20, help="Number of timesteps per mesh adaptation").tag(config=True)
    timestepper = Unicode('CrankNicolson').tag(config=True)

    # Solver
    params = PETScSolverParameters({'pc_type': 'lu',
                                    'mat_type': 'aij' ,
                                    'ksp_monitor': None,
                                    'ksp_converged_reason': None}).tag(config=True)

    # Physical 
    source_loc = List(default_value=None, allow_none=True).tag(config=True)
    source = FiredrakeScalarExpression(None, allow_none=True).tag(config=True)
    diffusivity = FiredrakeScalarExpression(Constant(1e-1)).tag(config=True)
    fluid_velocity = FiredrakeVectorExpression(None, allow_none=True).tag(config=True)

    # Adjoint
    solve_adjoint = Bool(False).tag(config=True)
    region_of_interest = List(default_value=[(20., 7.5, 0.5)]).tag(config=True)

    # Adaptivity
    h_min = PositiveFloat(1e-4, help="Minimum tolerated element size").tag(config=True)
    h_max = PositiveFloat(5., help="Maximum tolerated element size").tag(config=True)

    boundary_conditions = PETScSolverParameters({}).tag(config=True)

    def __init__(self, approach='fixed_mesh'):
        super(TracerOptions, self).__init__(approach)

        self.end_time -= 0.5*self.dt
        self.di = self.directory()  # TODO: surely this should be in __init__ anyway

    def set_diffusivity(self, fs):
        pass

    def set_velocity(self, fs):
        pass

    def set_source(self, fs):
        pass

    def set_initial_condition(self, fs):
        pass

    def set_kernel(self, fs):
        pass

    def exact_solution(self, fs):
        pass

    def exact_objective(self):
        return assemble(inner(self.kernel, self.solution)*dx)

class PowerOptions(TracerOptions):
    """
    Parameters for test case in [Power et al 2006].
    """
    def __init__(self, approach='fixed_mesh'):
        super(PowerOptions, self).__init__(approach)
        self.default_mesh = SquareMesh(40, 40, 4, 4)

        # Source / receiver
        self.source_loc = [(1., 2., 0.1)]
        self.region_of_interest = [(3., 2., 0.1)]

        # Boundary conditions  # TODO: make Thetis-conforming
        self.boundary_conditions[1] = 'dirichlet_zero'
        #self.boundary_conditions[2] = 'neumann_zero'  # FIXME
        self.boundary_conditions[3] = 'neumann_zero'
        self.boundary_conditions[4] = 'neumann_zero'

    def set_diffusivity(self, fs):
        self.diffusivity = Constant(1.)
        return self.diffusivity

    def set_velocity(self, fs):
        self.fluid_velocity = Function(fs)
        self.fluid_velocity.interpolate(as_vector((15., 0.)))
        return self.fluid_velocity

    def set_source(self, fs):
        self.source = Function(fs)
        #self.source.interpolate(self.bump(fs.mesh(), source=True))
        self.source.interpolate(self.box(fs.mesh(), source=True))
        area = assemble(self.source*dx)
        rescaling = 0.04/area if area != 0. else 1.
        self.source.interpolate(rescaling*self.source)
        self.source.rename("Source term")
        return self.source

    def set_objective_kernel(self, fs):
        self.kernel = Function(fs)
        #self.kernel.interpolate(self.bump(fs.mesh()))
        self.kernel.interpolate(self.box(fs.mesh()))
        area = assemble(self.kernel*dx)
        rescaling = 0.04/area if area != 0. else 1.
        self.kernel.interpolate(rescaling*self.kernel)
        return self.kernel


class TelemacOptions(TracerOptions):
    """
    Parameters for the 'Point source with diffusion' TELEMAC-2D test case.
    """
    def __init__(self, approach='fixed_mesh', offset=0.):
        super(TelemacOptions, self).__init__(approach)
        self.default_mesh = RectangleMesh(100, 20, 50, 10)
        self.offset = offset

        # Source / receiver
        self.source_loc = [(1.+self.offset, 5., 0.1)]
        self.region_of_interest = [(20., 7.5, 0.5)]
        self.source_value = 100.
        self.source_discharge = 0.1

        # Boundary conditions  # TODO: make Thetis-conforming
        self.boundary_conditions[1] = 'dirichlet_zero'
        self.boundary_conditions[3] = 'neumann_zero'
        self.boundary_conditions[4] = 'neumann_zero'

        # Time integration
        self.dt = 0.1
        self.end_time = 60.
        self.dt_per_export = 10
        self.dt_per_remesh = 20

    def set_diffusivity(self, fs):
        self.diffusivity = Constant(0.1)
        return self.diffusivity

    def set_velocity(self, fs):
        self.fluid_velocity = Function(fs)
        self.fluid_velocity.interpolate(as_vector((1., 0.)))
        return self.fluid_velocity

    def set_source(self, fs):
        x, y = SpatialCoordinate(fs.mesh())
        x0, y0, r0 = self.source_loc[0]
        self.source = Function(fs)
        with pyadjoint.stop_annotating():
            nrm=assemble(self.indicator(fs.mesh(), source=True)*dx)
        scaling = pi*r0*r0/nrm if nrm != 0 else 1
        scaling *= 0.5*self.source_value  # TODO: where does factor of half come from?
        self.source.interpolate(self.indicator(fs.mesh(), source=True, scale=scaling))
        return self.source

    def set_initial_condition(self, fs):
        self.initial_value = Function(fs)
        return self.initial_value

    def set_objective_kernel(self, fs):
        self.kernel = Function(fs)
        self.kernel.interpolate(self.indicator(fs.mesh()))
        return self.kernel

    def exact_solution(self, fs):
        self.solution = Function(fs)
        mesh = fs.mesh()
        x, y = SpatialCoordinate(mesh)
        x0, y0, r = self.source_loc[0]
        u = self.set_velocity(VectorFunctionSpace(fs.mesh(), fs.ufl_element()))
        nu = self.set_diffusivity(fs)
        #q = 0.01  # sediment discharge of source (kg/s)
        q = 1
        r = max_value(sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)), r)  # (Bessel fn explodes at (x0, y0))
        self.solution.interpolate(0.5*q/(pi*nu)*exp(0.5*u[0]*(x-x0)/nu)*bessk0(0.5*u[0]*r/nu))
        #self.solution.interpolate(0.5*q/(pi*nu)*exp(0.5*u[0]*(x-x0)/nu))
        #tmp = Function(fs).interpolate(0.5*u[0]*r/nu)
        #self.solution.dat.data[:] *= kn(0, tmp.dat.data)
        outfile = File(self.directory() + 'analytic.pvd')
        outfile.write(self.solution)  # NOTE: use 40 discretisation levels in ParaView
        return self.solution

    def exact_objective(self, fs1, fs2):
        self.exact_solution(fs1)
        self.set_objective_kernel(fs2)
        return assemble(self.kernel*self.solution*dx)


class TelemacOptions_Centred(TelemacOptions):
    def __init__(self, approach='fixed_mesh', offset=0.):
        super(TelemacOptions_Centred, self).__init__(approach, offset)
        self.region_of_interest = [(20., 5., 0.5)]


# TODO: Set three different tracers?
class LeVequeOptions(TracerOptions):
    """
    Parameters for test case in [LeVeque 1996].
    """
    def __init__(self, approach='fixed_mesh'):
        super(LeVequeOptions, self).__init__(approach)
        self.default_mesh = UnitSquareMesh(40, 40)

        # Source / receiver
        self.source_loc = [(0.25, 0.5, 0.15), (0.5, 0.25, 0.15), (0.5, 0.75, 0.15), (0.475, 0.525, 0.85)]
        self.region_of_interest = [(0.5, 0.75, 0.18)]

        # Boundary conditions
        q_in = Constant(1.0)
        for i in range(4):
            self.boundary_conditions[i] = {i: {'value': q_in}}

        # Time integration
        self.dt = math.pi/300.0
        self.end_time = 2*math.pi
        self.dt_per_export = 10
        self.dt_per_remesh = 10

        # Exact objective value
        #bell_r2 = self.source_loc[0][2]**2
        #cone_r2 = self.source_loc[1][2]**2
        cyl_x0, cyl_y0, cyl_r0 = self.source_loc[2]
        cyl_r2 = cyl_r0**2
        slot_left, slot_right, slot_top = self.source_loc[3]
        slot_left -= cyl_x0
        slot_right -= cyl_x0
        slot_top -= cyl_y0
        #bell = math.pi*bell_r2/4
        #cone = math.pi*cone_r2/3
        cyl = math.pi*cyl_r2
        slot = slot_top*(slot_right-slot_left)
        slot += -0.5*slot_right*math.sqrt(cyl_r2 - slot_right**2)
        slot += 0.5*slot_left*math.sqrt(cyl_r2 - slot_left**2)
        slot += -0.5*cyl_r2*math.atan(slot_right/math.sqrt(cyl_r2 - slot_right**2))
        slot += 0.5*cyl_r2*math.atan(slot_left/math.sqrt(cyl_r2 - slot_left**2))
        #self.J_exact = 1 + bell + cone + cyl - slot
        self.J_exact = 1 + cyl - slot
        #print("Exact objective: {:.4e}".format(self.J_exact))  # TODO: Check this

    def set_diffusivity(self, fs):
        self.diffusivity = Constant(0.)  # TODO: could use None?
        return self.diffusivity

    def set_velocity(self, fs):
        x, y = SpatialCoordinate(fs.mesh())
        self.fluid_velocity = Function(fs)
        self.fluid_velocity.interpolate(as_vector((0.5 - y, x - 0.5)))
        return self.fluid_velocity

    def set_initial_condition(self, fs):
        x, y = SpatialCoordinate(fs.mesh())
        bell_x0, bell_y0, bell_r0 = self.source_loc[0]
        cone_x0, cone_y0, cone_r0 = self.source_loc[1]
        cyl_x0, cyl_y0, cyl_r0 = self.source_loc[2]
        slot_left, slot_right, slot_top = self.source_loc[3]
        bell = 0.25*(1+cos(math.pi*min_value(sqrt(pow(x-bell_x0, 2) + pow(y-bell_y0, 2))/bell_r0, 1.0)))
        cone = 1.0 - min_value(sqrt(pow(x-cone_x0, 2) + pow(y-cone_y0, 2))/cyl_r0, 1.0)
        slot_cyl = conditional(sqrt(pow(x-cyl_x0, 2) + pow(y-cyl_y0, 2)) < cyl_r0,
                     conditional(And(And(x > slot_left, x < slot_right), y < slot_top),
                       0.0, 1.0), 0.0)

        self.initial_value = Function(fs)
        self.initial_value.interpolate(1.0 + bell + cone + slot_cyl)
        return self.initial_value

    def set_objective_kernel(self, fs):
        self.kernel = Function(fs)
        self.kernel.interpolate(self.indicator(fs.mesh()))
        area = assemble(self.kernel*dx)
        area_exact = math.pi*self.region_of_interest[0][2]**2
        rescaling = area_exact/area if area != 0. else 1
        self.kernel.interpolate(rescaling*self.kernel)
        return self.kernel
