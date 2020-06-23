from thetis import *
from thetis.physical_constants import *
from thetis.sediments import SedimentModel

from adapt_utils.solver import UnsteadyProblem
from adapt_utils.adapt.metric import *


__all__ = ["UnsteadyShallowWaterProblem"]


class UnsteadyShallowWaterProblem(UnsteadyProblem):
    """
    General solver object for time-dependent shallow water problems.
    """
    def __init__(self, op, mesh=None, **kwargs):
        p = op.degree
        if op.family == 'dg-dg' and p >= 0:
            fe = VectorElement("DG", triangle, p)*FiniteElement("DG", triangle, p)
        elif op.family == 'dg-cg' and p >= 0:
            fe = VectorElement("DG", triangle, p)*FiniteElement("Lagrange", triangle, p+1)
        else:
            raise NotImplementedError

        super(UnsteadyShallowWaterProblem, self).__init__(op, mesh, fe, **kwargs)
        prev_solution = kwargs.get('prev_solution')

        if prev_solution is not None:
            self.interpolate_solution(prev_solution)

        # Physical fields
        physical_constants['g_grav'].assign(op.g)

        # Classification
        self.nonlinear = True

    def set_fields(self):
        self.fields = {}
        self.fields['viscosity'] = self.op.set_viscosity(self.P1)
        self.fields['diffusivity'] = self.op.set_diffusivity(self.P1)
        self.fields['inflow'] = self.op.set_inflow(self.P1_vec)
        self.fields['coriolis'] = self.op.set_coriolis(self.P1)
        self.fields['quadratic_drag_coefficient'] = self.op.set_quadratic_drag_coefficient(self.P1DG)
        self.fields['manning_drag_coefficient'] = self.op.set_manning_drag_coefficient(self.P1)
        self.fields['source'] = self.op.source

        self.op.set_boundary_surface()

    def create_solutions(self):
        super(UnsteadyShallowWaterProblem, self).create_solutions()
        u, eta = self.solution.split()
        u.rename("Fluid velocity")
        eta.rename("Elevation")
        z, zeta = self.adjoint_solution.split()
        z.rename("Adjoint fluid velocity")
        zeta.rename("Adjoint elevation")

    def set_stabilisation(self):
        self.stabilisation = self.stabilisation or 'no'
        if self.stabilisation in ('no', 'lax_friedrichs'):
            self.stabilisation_parameter = self.op.stabilisation_parameter
        else:
            raise ValueError("Stabilisation method {:s} for {:s} not recognised".format(self.stabilisation, self.__class__.__name__))

    def get_tracer(self, adjoint=False):
        """
        Retrieve forward or adjoint solution, as specified by boolean kwarg `adjoint`.
        """
        return self.adjoint_tracer if adjoint else self.solution_old_tracer

    def get_bathymetry(self, adjoint=False):
        """
        Retrieve forward or adjoint solution, as specified by boolean kwarg `adjoint`.
        """
        return self.adjoint_bath if adjoint else self.solution_old_bathymetry

    def project_tracer(self, val, adjoint=False):
        """
        Project forward or adjoint solution, as specified by the boolean kwarg
        `adjoint`.
        """
        self.project(val, out=self.get_tracer(adjoint=adjoint))

    def project_bathymetry(self, val, adjoint=False):
        """
        Project forward or adjoint solution, as specified by the boolean kwarg
        `adjoint`.
        """
        self.project(val, out=self.get_bathymetry(adjoint=adjoint))

    def solve_step(self, adjoint=False):
        try:
            assert not adjoint
        except AssertionError:
            raise NotImplementedError  # TODO
        if not hasattr(self, 'solver_obj'):
            self.setup_solver_forward()

        if self.op.suspended:
            if self.solution_old_tracer is not None:
                self.tracer_interp = project(self.solution_old_tracer, self.P1DG)
            else:
                self.tracer_interp = project(self.solver_obj.fields.sediment_2d, self.P1DG)
                self.solution_old_tracer = project(self.solver_obj.fields.sediment_2d, self.P1DG)

        u_interp, eta_interp = self.solution.split()

        if self.op.suspended:
            if self.op.tracer_init is not None:
                self.solver_obj.assign_initial_conditions(uv=u_interp, elev=eta_interp, sediment=self.tracer_interp)

        self.solver_obj.fields.bathymetry_2d.project(self.solution_old_bathymetry)

        self.solver_obj.options.simulation_end_time = self.step_end - 0.5*self.op.dt
        #self.op.update(self.solver_obj.options)

        self.solver_obj.iterate(update_forcings=self.op.get_update_forcings(self.solver_obj),
                                export_func=self.op.get_export_func(self.solver_obj))

        old_mesh = Mesh(Function(self.mesh.coordinates))
        P1DG_old = FunctionSpace(old_mesh, "DG", 1)
        P1_old = FunctionSpace(old_mesh, "CG", 1)
        
        self.solution = self.solver_obj.fields.solution_2d

        if isinstance(self.solver_obj.fields.bathymetry_2d, Constant):
            solution_bathymetry = Constant(self.solver_obj.fields.bathymetry_2d)
            self.solution_old_bathymetry = Constant(solution_bathymetry)
        else:
            solution_bathymetry = self.solver_obj.fields.bathymetry_2d.copy(deepcopy=True)
            self.solution_old_bathymetry = project(solution_bathymetry, P1_old)

        if self.op.suspended:
            solution_tracer = self.solver_obj.fields.sediment_2d.copy(deepcopy=True)
            self.solution_old_tracer = project(solution_tracer, P1DG_old)

    def setup_solver_forward(self):
        if not hasattr(self, 'remesh_step'):
            self.remesh_step = 0
        op = self.op

        # Use appropriate bathymetry
        if hasattr(self, 'solution_old_bathymetry') and self.solution_old_bathymetry is not None:
            if isinstance(self.solution_old_bathymetry, Constant):
                op.bathymetry = Constant(self.solution_old_bathymetry)
            else:
                op.bathymetry = project(self.solution_old_bathymetry, self.P1)
        else:
            op.bathymetry = self.op.set_bathymetry(self.P1)
            self.solution_old_bathymetry = self.op.set_bathymetry(self.P1)
        b = op.bathymetry
        self.solver_obj = solver2d.FlowSolver2d(self.mesh, b)

        self.solver_obj.export_initial_state = self.remesh_step == 0

        # Initial conditions
        u_interp, eta_interp = self.solution.split()
        if op.suspended:
            if hasattr(self, 'solution_old_tracer') and self.solution_old_tracer is not None:
                self.tracer_interp = project(self.solution_old_tracer, self.P1DG)
            else:
                self.tracer_interp = project(self.op.tracer_init, self.P1DG)
                self.solution_old_tracer = project(self.op.tracer_init, self.P1DG)

            self.uv_d, self.eta_d = self.solution.split()    

        options = self.solver_obj.options

        self.solver_obj.sediment_model = op.sed_mod
        
        options.update(self.solver_obj.sediment_model.options)
        
        options.use_nonlinear_equations = self.nonlinear
        options.check_volume_conservation_2d = True
        if hasattr(options, 'use_lagrangian_formulation'):  # TODO: Temporary
            options.use_lagrangian_formulation = op.approach == 'ale'

        # Timestepping
        options.timestep = op.dt
        options.simulation_export_time = op.dt*op.dt_per_export
        options.simulation_end_time = self.step_end - 0.5*op.dt
        options.timestepper_type = op.timestepper
        if op.params != {}:
            options.timestepper_options.solver_parameters = op.params
        if op.debug:
            # options.timestepper_options.solver_parameters['snes_monitor'] = None
            print_output(options.timestepper_options.solver_parameters)
        if hasattr(options.timestepper_options, 'implicitness_theta'):
            options.timestepper_options.implicitness_theta = op.implicitness_theta

        # Outputs
        options.output_directory = self.di

        if op.suspended:
            options.fields_to_export = ['uv_2d', 'elev_2d', 'sediment_2d'] if op.plot_pvd else []
            options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d', 'sediment_2d'] if op.save_hdf5 else []
        else:
            options.fields_to_export = ['uv_2d', 'elev_2d'] if op.plot_pvd else []
            options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d'] if op.save_hdf5 else []

        # Parameters
        options.use_grad_div_viscosity_term = op.grad_div_viscosity
        options.element_family = op.family
        options.horizontal_viscosity = self.fields['viscosity']
        options.horizontal_diffusivity = self.fields['diffusivity']
        if op.friction == 'nikuradse':
            options.nikuradse_bed_roughness = op.sed_mod.ksp
        elif op.friction == 'manning':
            options.manning_drag_coefficient = self.fields['manning_drag_coefficient']
        options.coriolis_frequency = self.fields['coriolis']
        options.use_lax_friedrichs_velocity = self.stabilisation == 'lax_friedrichs'
        options.lax_friedrichs_velocity_scaling_factor = self.stabilisation_parameter
        options.morphological_acceleration_factor = op.morfac
        if hasattr(op, 'bnd_dict'):
            options.equilibrium_sediment_bd_ids = op.bnd_dict
        options.use_grad_depth_viscosity_term = op.grad_depth_viscosity
        options.use_automatic_sipg_parameter = op.sipg_parameter is None
        options.use_wetting_and_drying = op.wetting_and_drying
        options.wetting_and_drying_alpha = op.wetting_and_drying_alpha

        if hasattr(op, 'norm_smoother_constant'):
            options.norm_smoother = op.norm_smoother_constant      

        if op.suspended:
            if op.depth_integrated and op.conservative:
                options.tracer_depth_integ_source = op.sed_mod.ero
                options.tracer_depth_integ_sink = op.sed_mod.depo_term
            else:
                options.tracer_source_2d = op.sed_mod.ero_term
                options.tracer_sink_2d = op.sed_mod.depo_term

        options.solve_exner = True
        # Boundary conditions
        self.solver_obj.bnd_functions['shallow_water'] = op.set_boundary_conditions(self.V)

        if op.suspended:
            self.solver_obj.bnd_functions['sediment'] = op.set_boundary_conditions_tracer(op.sed_mod)

            for i in self.solver_obj.bnd_functions['sediment'].keys():
                if i in self.solver_obj.bnd_functions['shallow_water'].keys():
                    self.solver_obj.bnd_functions['sediment'][i].update(self.solver_obj.bnd_functions['shallow_water'][i])
            for i in self.solver_obj.bnd_functions['shallow_water'].keys():
                if i not in self.solver_obj.bnd_functions['sediment'].keys():
                    self.solver_obj.bnd_functions['sediment'].update({i:self.solver_obj.bnd_functions['shallow_water'][i]})        

        if op.suspended:
            if self.op.tracer_init is not None:
                self.solver_obj.assign_initial_conditions(uv=u_interp, elev=eta_interp, sediment=self.tracer_interp)
        else:
            self.solver_obj.assign_initial_conditions(uv=u_interp, elev=eta_interp)

        if hasattr(self, 'extra_setup'):
            self.extra_setup()
        # Ensure correct iteration count
        self.solver_obj.i_export = self.remesh_step
        self.solver_obj.next_export_t = self.remesh_step*op.dt*op.dt_per_remesh
        self.solver_obj.iteration = int(self.remesh_step*op.dt_per_remesh)
        self.solver_obj.simulation_time = self.remesh_step*op.dt*op.dt_per_remesh
        for e in self.solver_obj.exporters.values():
            e.set_next_export_ix(self.solver_obj.i_export)

    def get_bnd_functions(self, *args):
        b = self.op.bathymetry if options.solve_sediment else self.fields['bathymetry']
        swt = shallowwater_eq.ShallowWaterTerm(self.V, bathymetry=b)
        return swt.get_bnd_functions(self.boundary_conditions, *args)

    def get_qoi_kernel(self):
        self.kernel = self.op.set_qoi_kernel(self.solver_obj)

    def plot_solution(self, adjoint=False):
        if adjoint:
            z, zeta = self.adjoint_solution.split()
            z.rename("Adjoint fluid velocity")
            zeta.rename("Adjoint elevation")
            self.adjoint_solution_file.write(z, zeta)
        else:
            u, eta = self.solution.split()
            u.rename("Fluid velocity")
            eta.rename("Elevation")
            self.solution_file.write(u, eta)

    def get_hessian_metric(self, adjoint=False, **kwargs):
        kwargs.setdefault('noscale', False)
        kwargs.setdefault('degree', 1)
        kwargs.setdefault('mesh', self.mesh)
        kwargs.setdefault('op', self.op)
        field = self.op.adapt_field
        sol = self.adjoint_solution if adjoint else self.solution
        u, eta = sol.split()

        fdict = {'elevation': eta, 'velocity_x': u[0], 'velocity_y': u[1],
                 'speed': sqrt(inner(u, u)), 'inflow': inner(u, self.fields['inflow'])}
        fdict.update(self.fields)

        self.M = Function(self.P1_ten)
        if field in fdict:
            self.M = steady_metric(fdict[field], **kwargs)
        elif field == 'all_avg':
            self.M += steady_metric(fdict['velocity_x'], **kwargs)
            self.M += steady_metric(fdict['velocity_y'], **kwargs)
            self.M += steady_metric(fdict['elevation'], **kwargs)
            self.M /= 3.0
        elif field == 'all_int':
            self.M = metric_intersection(steady_metric(fdict['velocity_x'], **kwargs),
                                         steady_metric(fdict['velocity_y'], **kwargs))
            self.M = metric_intersection(self.M, steady_metric(fdict['elevation'], **kwargs))
        elif 'avg' in field and 'int' in field:
            raise NotImplementedError  # TODO
        elif 'avg' in field:
            fields = field.split('_avg_')
            num_fields = len(fields)
            for i in range(num_fields):
                self.M += steady_metric(fdict[fields[i]], **kwargs)
            self.M /= num_fields
        elif 'int' in field:
            fields = field.split('_int_')
            self.M = steady_metric(fdict[fields[0]], **kwargs)
            for i in range(1, len(fields)):
                self.M = metric_intersection(self.M, steady_metric(fields[f[i]], **kwargs))
        else:
            raise ValueError("Adaptation field {:s} not recognised.".format(field))

    def custom_adapt(self):
        if self.approach == 'vorticity':
            self.indicator = Function(self.P1, name='vorticity')
            self.indicator.interpolate(curl(self.solution.split()[0]))
            self.get_isotropic_metric()
