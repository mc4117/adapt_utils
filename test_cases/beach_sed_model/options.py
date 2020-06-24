from thetis import *
from thetis.configuration import *
from thetis.options import ModelOptions2d
from thetis.sediments import SedimentModel

from adapt_utils.swe.morphological.morphological_options import MorphOptions

import os
import time
import datetime
import numpy as np
from matplotlib import rc

rc('text', usetex=True)


__all__ = ["BeachOptions"]


class BeachOptions(MorphOptions):
    """
    Parameters for test case adapted from [1].

    [1] Roberts, W. et al. "Investigation using simple mathematical models of 
    the effect of tidal currents and waves on the profile shape of intertidal 
    mudflats." Continental Shelf Research 20.10-11 (2000): 1079-1097.
    """

    def __init__(self, friction='manning', plot_timeseries=False, nx=1, ny=1, mesh = None, input_dir = None, output_dir = None, **kwargs):

        super(BeachOptions, self).__init__(**kwargs)
        
        try:
            assert friction in ('nikuradse', 'manning')
        except AssertionError:
            raise ValueError("Friction parametrisation '{:s}' not recognised.".format(friction))
        self.friction = friction  
        
        self.lx = 220
        self.ly = 10

        if output_dir is not None:
            self.di = output_dir
        else:
            self.export_intermediate = False

        self.plot_timeseries = plot_timeseries
        if mesh is None:
            self.default_mesh = RectangleMesh(np.int(220*nx), 10*ny, self.lx, self.ly)
            self.P1DG = FunctionSpace(self.default_mesh, "DG", 1)
            self.P1 = FunctionSpace(self.default_mesh, "CG", 1)
            self.P1_vec = VectorFunctionSpace(self.default_mesh, "CG", 1)
            self.P1_vec_dg = VectorFunctionSpace(self.default_mesh, "DG", 1)
        else:
            self.P1DG = FunctionSpace(mesh, "DG", 1)
            self.P1 = FunctionSpace(mesh, "CG", 1)
            self.P1_vec = VectorFunctionSpace(mesh, "CG", 1)
            self.P1_vec_dg = VectorFunctionSpace(mesh, "DG", 1)
            
        if input_dir is not None:
            self.input_dir = input_dir

        self.plot_pvd = True
        self.hessian_recovery = 'dL2'

        self.grad_depth_viscosity = True

        if self.export_intermediate:
            self.bathymetry_file = File(self.di + "/bathy.pvd")

        self.num_hours = 3

        self.t_old = Constant(0.0)
        # Stabilisation
        self.stabilisation = 'lax_friedrichs'
        
        self.morfac = Constant(25)

        if mesh is None:
            self.set_up_morph_model(self.default_mesh)
        else:
            self.set_up_morph_model(mesh)

        # Boundary conditions
        h_amp = 0.5  # Ocean boundary forcing amplitude
        v_amp = 1 # Ocean boundary foring velocity
        omega = 0.5  # Ocean boundary forcing frequency
        self.ocean_elev_func = lambda t: (h_amp * np.cos(-omega *(t+(100.0))))
        self.ocean_vel_func = lambda t: (v_amp * np.cos(-omega *(t+(100.0))))
        
        self.tracer_init = Constant(0.0)    

        # Time integration

        self.dt = 0.05
        self.end_time = float(self.num_hours*3600.0/self.morfac)
        self.dt_per_export = 80
        self.dt_per_remesh = 80
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

        # Outputs
        if self.export_intermediate:
            self.bath_file = File(os.path.join(self.di, 'bath_export.pvd'))
            self.diff_bath_file = File(os.path.join(self.di, 'diff_export.pvd'))
            #self.eta_tilde_file = File(os.path.join(self.di, 'eta_tilde.pvd'))
        self.eta_tilde = Function(self.P1DG, name='Modified elevation')        

        # Physical
        self.base_viscosity = 1 #0.2 #5*10**(-1)
        self.base_diffusivity = 0.15
        self.gravity = Constant(9.81)
        self.porosity = Constant(0.4)
        self.ks = Constant(0.025)
        self.average_size = 0.0002  # Average sediment size
        self.ksp = Constant(3*self.average_size)

        self.wetting_and_drying = True
        self.depth_integrated = True
        self.conservative = True
        self.slope_eff = True
        self.angle_correction = False
        self.convective_vel_flag = True
        self.solve_tracer = True 

        self.suspended = True
        self.bedload = True

        self.wetting_and_drying_alpha = Constant(5/25)
        self.norm_smoother_constant = Constant(10/25)

        # Initial

        self.elev_init, self.uv_init = self.initialise_fields(self.input_dir, self.di)
        
        self.uv_d = Function(self.P1_vec_dg).project(self.uv_init)

        self.eta_d = Function(self.P1DG).project(self.elev_init)
        
        if not hasattr(self, 'bathymetry') or self.bathymetry is None:
            self.bathymetry = self.set_bathymetry(self.P1)
            self.bath_init = self.set_bathymetry(self.P1)
            self.diff_bathy = Function(self.P1)

        self.sed_mod = SedimentModel(ModelOptions2d(), suspendedload=self.suspended, convectivevel=self.convective_vel_flag,
                            bedload=self.bedload, angle_correction=self.angle_correction, slope_eff=self.slope_eff, seccurrent=False,
                            mesh2d=mesh, bathymetry_2d=self.bathymetry,
                            uv_init = self.uv_d, elev_init = self.eta_d, ks=self.ks, average_size=self.average_size, 
                            cons_tracer = self.conservative, wetting_and_drying = self.wetting_and_drying, wetting_alpha = self.wetting_and_drying_alpha)


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
        self.bathymetry.interpolate(4 - x/25)
        return self.bathymetry

    def set_viscosity(self, fs):
        x, y = SpatialCoordinate(fs.mesh())
        self.viscosity = Function(fs)
        sponge_viscosity = Function(fs).interpolate(conditional(x>=100, -99 + x, Constant(1.0)))
        self.viscosity.interpolate(sponge_viscosity*self.base_viscosity)
        return self.viscosity

    def set_coriolis(self, fs):
        return

    def set_boundary_conditions(self, fs):
        if not hasattr(self, 'elev_in'):
            self.elev_in = Constant(0.0)
        if not hasattr(self, 'vel_in'):
            self.vel_in = Constant(as_vector((0.0, 0.0)))
        self.elev_in.assign(self.ocean_elev_func(0.0))
        vel_const = Constant(self.ocean_vel_func(0.0))
        self.vel_in.assign(as_vector((vel_const, 0.0)))

        inflow_tag = 1
        boundary_conditions = {}
        boundary_conditions[inflow_tag] = {'elev': self.elev_in, 'uv': self.vel_in}
        return boundary_conditions

    def update_boundary_conditions(self, solver_obj, t=0.0):
        self.elev_in.assign(self.ocean_elev_func(t))
        vel_const = Constant(self.ocean_vel_func(t))
        self.vel_in.assign(as_vector((vel_const, 0.0)))

    def set_initial_condition(self, fs):
        """
        Set initial elevation and velocity using asymptotic solution.

        :arg fs: `FunctionSpace` in which the initial condition should live.
        """
        self.initial_value = Function(fs, name="Initial condition")
        u, eta = self.initial_value.split()
        u.project(self.uv_init)
        eta.project(self.elev_init)
        return self.initial_value

    def set_boundary_conditions_tracer(self, fs):
        inflow_tag = 1
        boundary_conditions = {}
        #boundary_conditions[inflow_tag] = {'value': self.sed_mod.equiltracer}
        return boundary_conditions

    def get_update_forcings(self, solver_obj):

        def update_forcings(t):

            self.update_boundary_conditions(solver_obj, t=t)

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
            #File(outputdir + "/elevation_imported.pvd").write(elev_init)
            chk.close()
        # velocity
        with timed_stage('initialising velocity'):
            chk = DumbCheckpoint(inputdir + "/velocity" , mode=FILE_READ)
            uv_init = Function(vector_dg, name="velocity")
            chk.load(uv_init)
            #File(outputdir + "/velocity_imported.pvd").write(uv_init)
            chk.close()

        return  elev_init, uv_init, 

    def get_export_func(self, solver_obj):
        self.bath_export = solver_obj.fields.bathymetry_2d
        bathymetry_displacement = solver_obj.depth.wd_bathymetry_displacement
        eta = solver_obj.fields.elev_2d        

        def export_func():
            self.eta_tilde.project(eta + bathymetry_displacement(eta))
            #self.eta_tilde_file.write(self.eta_tilde)
            self.bath_file.write(self.bath_export)
            self.diff_bathy.project(self.bath_export - self.bath_init)
            self.diff_bath_file.write(self.diff_bathy)
        return export_func 

    def set_boundary_surface(self):
        """Set the initial displacement of the boundary elevation."""
        pass
