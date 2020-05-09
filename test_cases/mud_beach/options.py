from thetis import *
from thetis.configuration import *

from adapt_utils.swe.morphological_options import MorphOptions

import os
import time
import datetime
import numpy as np
from matplotlib import rc

rc('text', usetex=True)


__all__ = ["MudBeachOptions"]


class MudBeachOptions(MorphOptions):
    """
    Parameters for test case adapted from [1].

    [1] Roberts, W. et al. "Investigation using simple mathematical models of 
    the effect of tidal currents and waves on the profile shape of intertidal 
    mudflats." Continental Shelf Research 20.10-11 (2000): 1079-1097.
    """

    def __init__(self, friction='manning', plot_timeseries=False, nx=1, ny=1, mesh = None, **kwargs):
        
        self.settling_velocity = Constant(0.0005)
        super(MudBeachOptions, self).__init__(**kwargs)
        
        try:
            assert friction in ('nikuradse', 'manning')
        except AssertionError:
            raise ValueError("Friction parametrisation '{:s}' not recognised.".format(friction))
        self.friction = friction     
        
        self.lx = 10000
        self.ly = 800

        self.plot_timeseries = plot_timeseries
        if mesh is None:
            self.default_mesh = RectangleMesh(np.int(50*nx), 4*ny, self.lx, self.ly)
            self.P1DG = FunctionSpace(self.default_mesh, "DG", 1)
            self.P1 = FunctionSpace(self.default_mesh, "CG", 1)
            self.P1_vec = VectorFunctionSpace(self.default_mesh, "CG", 1)
            self.P1_vec_dg = VectorFunctionSpace(self.default_mesh, "DG", 1)
        else:
            print('passed here')
            self.P1DG = FunctionSpace(mesh, "DG", 1)
            self.P1 = FunctionSpace(mesh, "CG", 1)
            self.P1_vec = VectorFunctionSpace(mesh, "CG", 1)
            self.P1_vec_dg = VectorFunctionSpace(mesh, "DG", 1)

        self.plot_pvd = True
        self.hessian_recovery = 'dL2'

        self.grad_depth_viscosity = True

        self.bathymetry_file = File(self.di + "/bathy.pvd")

        self.num_hours = 3650*24

        self.t_old = Constant(0.0)
        # Stabilisation
        self.stabilisation = 'lax_friedrichs'
        
        self.morfac = 365

        if mesh is None:
            self.set_up_morph_model()
        else:
            self.set_up_morph_model(mesh)

        # Boundary conditions
        h_amp = -2  # Ocean boundary forcing amplitude
        v_amp = -0.4 # Ocean boundary foring velocity
        h_T = 12*3600.  # Ocean boundary forcing period
        self.ocean_elev_func = lambda t: h_amp * sin(2 * pi * t / h_T)
        self.ocean_vel_func = lambda t: v_amp * cos(2 * pi * t / h_T)
        self.flux_func = lambda t, d: -self.ocean_vel_func(t)*d.at([0.0, self.ly/2])*self.ly
        
        self.tracer_func = lambda t: -((10**-4)/2.65)*cos(2*pi*t/(12*3600)) if 10800 <= t%(12*3600) < 32400 else 0.0

        # Time integration

        self.dt = 75
        self.end_time = self.num_hours*3600.0/self.morfac
        self.dt_per_export = 200
        self.dt_per_remesh = 200
        self.timestepper = 'CrankNicolson'
        self.implicitness_theta = 1.0       
        
        # Adaptivity
        self.h_min = 1e-8
        self.h_max = 10.

        # Goal-Oriented
        self.qoi_mode = 'inundation_volume'

        # Timeseries
        self.wd_obs = []
        self.trange = np.linspace(0.0, self.end_time, self.num_hours+1)
        tol = 1e-8  # FIXME: Point evaluation hack
        self.xrange = np.linspace(tol, 16-tol, 20)
        self.qois = []

        
    def set_up_morph_model(self, mesh = None):

        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        outputdir = 'outputs' + st

        self.di = outputdir  # "morph_output"
        # Outputs
        self.bath_file = File(os.path.join(self.di, 'bath_export.pvd'))
        self.eta_tilde_file = File(os.path.join(self.di, 'eta_tilde.pvd'))
        self.eta_tilde = Function(self.P1DG, name='Modified elevation')        

        # Physical
        self.base_viscosity = 1e-6        
        self.base_diffusivity = 0.15
        self.gravity = Constant(9.81)
        self.porosity = Constant(0.4)
        self.ks = Constant(0.025)
        self.average_size = 150*(10**(-6))  # Average sediment size
        self.ksp = Constant(3*self.average_size)

        self.wetting_and_drying = True
        self.depth_integrated = True
        self.conservative = True
        self.implicit_source = True
        self.solve_tracer = True 
        self.slope_eff = True
        self.angle_correction = True
        
        self.wetting_and_drying_alpha = Constant(0.2)
        self.norm_smoother_constant = Constant(0.2)

        # Initial
        self.elev_init = Function(self.P1DG).interpolate(Constant(0.0))
        self.uv_init = as_vector((10**(-7), 0.0))
        
        self.uv_d = Function(self.P1_vec_dg).interpolate(self.uv_init)

        self.eta_d = Function(self.P1DG).interpolate(self.elev_init)
        
        self.fixed_tracer = Constant(0.0)
        
        if mesh is None:
            self.set_up_suspended(self.default_mesh, tracer = self.fixed_tracer)
            self.set_up_bedload(self.default_mesh)
        else:
            self.set_up_suspended(mesh, tracer = self.fixed_tracer)
            self.set_up_bedload(mesh)
        
    def set_source_tracer(self, fs, solver_obj=None, init=False):
        if init:
            self.depo = Function(fs).interpolate(self.settling_velocity * self.coeff)
            self.ero = Function(fs).interpolate(self.settling_velocity * self.ceq)
        else:
            self.depo.interpolate(self.settling_velocity * self.coeff)
            self.ero.interpolate(self.settling_velocity * self.ceq)
        return self.depo, self.ero

    def set_quadratic_drag_coefficient(self, fs):
        if self.friction == 'nikuradse':
            self.quadratic_drag_coefficient = project(self.get_cfactor(), fs)
        return self.quadratic_drag_coefficient

    def get_cfactor(self):
        return Function(self.P1DG).interpolate(Constant(0.002))

    def set_manning_drag_coefficient(self, fs):
        if self.friction == 'manning':
            if hasattr(self, 'friction_coeff'):
                self.manning_drag_coefficient = Constant(self.friction_coeff)
            else:
                self.manning_drag_coefficient = Constant(0.02)
        return self.manning_drag_coefficient

    def set_bathymetry(self, fs, **kwargs):

        x, y = SpatialCoordinate(fs.mesh())
        self.bathymetry = Function(fs, name="Bathymetry")
        self.bathymetry.interpolate(4.5 - x/1000)
        return self.bathymetry

    def set_viscosity(self, fs):
        self.viscosity = Function(fs)
        self.viscosity.assign(self.base_viscosity)
        return self.viscosity

    def set_coriolis(self, fs):
        return

    def set_boundary_conditions(self, fs):
        if not hasattr(self, 'flux_in'):
            self.set_boundary_flux()
        if not hasattr(self, 'depth'):
            d_init = Function(self.elev_init.function_space).interpolate(self.elev_init + self.bathymetry.at([0,0]))
        self.flux_in.assign(self.flux_func(0.0, self.depth))
        inflow_tag = 1
        boundary_conditions = {}
        boundary_conditions[inflow_tag] = {'flux': self.flux_in,}
        return boundary_conditions

    def update_boundary_conditions(self, solver_obj, t=0.0):
        depth = Function(self.P1DG).interpolate(solver_obj.depth.get_total_depth(self.ocean_elev_func(t)))
        self.flux_in.assign(self.flux_func(t, depth))

    def set_initial_condition(self, fs):
        self.initial_value = Function(fs, name="Initial condition")
        u, eta = self.initial_value.split()
        u.interpolate(as_vector([1.0e-7, 0.0]))
        eta.assign(0.0)
        return self.initial_value

    def set_boundary_conditions_tracer(self, fs):
        if not hasattr(self, 'tracer_flux_value'):
            self.tracer_flux_value = Constant(0.0)
        self.tracer_flux_value.assign(self.tracer_func(0.0))
        inflow_tag = 1
        boundary_conditions = {}
        boundary_conditions[inflow_tag] = {'value': self.tracer_flux_value}
        return boundary_conditions
    
    def update_boundary_conditions_tracer(self, t=0.0):
        self.tracer_flux_value.assign(self.tracer_func(t))

    def get_update_forcings(self, solver_obj):

        def update_forcings(t):
            
            self.update_boundary_conditions(solver_obj, t=t)
            self.update_boundary_conditions_tracer(t=t)
          
            self.update_key_hydro(solver_obj)
            
            if self.t_old.dat.data[:] == t:
                if self.suspended:
                    self.update_suspended(solver_obj)
                if self.bedload:
                    self.update_bedload(solver_obj)

                solve(self.f == 0, self.z_n1)

                self.bathymetry.assign(self.z_n1)
                solver_obj.fields.bathymetry_2d.assign(self.z_n1)
            
            self.t_old.assign(t)

        return update_forcings

    def initialise_fields(self, inputdir, outputdir):
        """
        Initialise simulation with results from a previous simulation
        """     
        from firedrake.petsc import PETSc
        try:
            import firedrake.cython.dmplex as dmplex
        except:
            import firedrake.dmplex as dmplex  # Older version        
        # mesh
        with timed_stage('mesh'):
            # Load
            newplex = PETSc.DMPlex().create()
            newplex.createFromFile(inputdir + '/myplex.h5')
            mesh = Mesh(newplex)
    
        DG_2d = FunctionSpace(mesh, 'DG', 1)  
        vector_dg = VectorFunctionSpace(mesh, 'DG', 1)          
        # elevation
        with timed_stage('initialising elevation'):
            chk = DumbCheckpoint(inputdir + "/elevation", mode=FILE_READ)
            elev_init = Function(DG_2d, name="elevation")
            chk.load(elev_init)
            File(outputdir + "/elevation_imported.pvd").write(elev_init)
            chk.close()
        # velocity
        with timed_stage('initialising velocity'):
            chk = DumbCheckpoint(inputdir + "/velocity" , mode=FILE_READ)
            uv_init = Function(vector_dg, name="velocity")
            chk.load(uv_init)
            File(outputdir + "/velocity_imported.pvd").write(uv_init)
            chk.close()

        return  elev_init, uv_init, 

    def get_export_func(self, solver_obj):
        self.bath_export = solver_obj.fields.bathymetry_2d
        bathymetry_displacement = solver_obj.depth.wd_bathymetry_displacement
        eta = solver_obj.fields.elev_2d        

        def export_func():
            self.eta_tilde.project(eta + bathymetry_displacement(eta))
            self.eta_tilde_file.write(self.eta_tilde) 
            self.bath_file.write(self.bath_export)
        return export_func 

    def set_boundary_surface(self):
        """Set the initial displacement of the boundary elevation."""
        pass
