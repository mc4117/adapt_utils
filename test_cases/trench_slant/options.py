from thetis import *
from thetis.configuration import *

from adapt_utils.swe.morphological.morphological_options import MorphOptions
from thetis.sediments import SedimentModel
from thetis.options import ModelOptions2d

import os
import time
import datetime
import numpy as np
from matplotlib import rc

rc('text', usetex=True)


__all__ = ["TrenchSlantOptions"]


class TrenchSlantOptions(MorphOptions):
    """
    Parameters for test case described in [1].

    [1] Clare, Mariana, et al. “Hydro-morphodynamics 2D Modelling Using a Discontinuous Galerkin Discretisation.”
    EarthArXiv, 9 Jan. 2020. Web.
    """

    def __init__(self, friction='manning', plot_timeseries=False, nx=1, ny=1, mesh = None, input_dir = None, output_dir = None, **kwargs):
        
        super(TrenchSlantOptions, self).__init__(**kwargs)
        
        try:
            assert friction in ('nikuradse', 'manning')
        except AssertionError:
            raise ValueError("Friction parametrisation '{:s}' not recognised.".format(friction))
        self.friction = friction        


        self.plot_timeseries = plot_timeseries
        if mesh is None:
            self.t_old = Constant(0.0)
            self.default_mesh = RectangleMesh(np.int(16*5*nx), np.int(np.ceil(5*ny)), 16, 1.1)
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
            
        if output_dir is not None:
            self.di = output_dir
        else:
            self.export_intermediate = False            

        self.plot_pvd = True
        self.hessian_recovery = 'dL2'

        self.grad_depth_viscosity = True

        self.num_hours = 15

        self.t_old = Constant(0.0)
        # Stabilisation
        self.stabilisation = 'lax_friedrichs'
        
        self.morfac = Constant(100)
        self.norm_smoother_constant = Constant(0.1)

        if mesh is None:
            self.set_up_morph_model(self.input_dir, self.default_mesh)
        else:
            self.set_up_morph_model(self.input_dir, mesh)

        # Time integration
        self.dt = 0.2
        self.end_time = float(self.num_hours*3600.0/self.morfac)
        self.dt_per_export = 40
        self.dt_per_remesh = 40
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
        
        #self.bnd_dict = {1}


        
    def set_up_morph_model(self, input_dir, mesh = None):

        # Outputs
        if self.export_intermediate:
            self.bath_file = File(os.path.join(self.di, 'bath_export.pvd'))     

        # Physical
        self.base_viscosity = 1e-6        
        self.base_diffusivity = 0.18161630470135287

        self.porosity = Constant(0.4)
        self.ks = Constant(0.025)
        self.average_size = 160*(10**(-6))  # Average sediment size        

        self.wetting_and_drying = False
        self.conservative = False
        self.implicit_source = False
        self.slope_eff = True
        self.angle_correction = True
        self.solve_tracer = True
        self.suspended = True
        self.convectivevel_flag = True
        self.bedload = True

        # Initial
        self.elev_init, self.uv_init = self.initialise_fields(mesh, input_dir, self.di)

        #self.uv_d = Function(self.P1_vec_dg).project(self.uv_init)

        #self.eta_d = Function(self.P1DG).project(self.elev_init)

        if not hasattr(self, 'bathymetry') or self.bathymetry is None:
            self.bathymetry = self.set_bathymetry(self.P1)


        if self.suspended:
            self.tracer_init = None


    def set_bathymetry(self, fs, **kwargs):

        initialdepth = Constant(0.297)
        depth_riv = Constant(initialdepth - 0.397)
        depth_trench = Constant(depth_riv - 0.15)
        depth_diff = depth_trench - depth_riv
        x, y = SpatialCoordinate(fs.mesh())
        trench = conditional(le(x, 5), (0.1*(y-0.55)) + depth_riv, conditional(le(x, 6.5), (0.1*(y-0.55)) + (1/1.5)*depth_diff*(x-6.5) + depth_trench,
                        conditional(le(x, 9.5), (0.1*(y-0.55)) + depth_trench, conditional(le(x, 11), (0.1*(y-0.55)) - (1/1.5)*depth_diff*(x-11) + depth_riv, (0.1*(y-0.55)) + depth_riv))))
        self.bathymetry = Function(fs, name="Bathymetry")
        self.bathymetry.interpolate(-trench)
        return self.bathymetry

    def set_viscosity(self, fs):
        self.viscosity = Function(fs)
        self.viscosity.assign(self.base_viscosity)
        return self.viscosity

    def set_coriolis(self, fs):
        return

    def set_boundary_conditions(self, fs):
        inflow_tag = 1
        outflow_tag = 2
        boundary_conditions = {}
        boundary_conditions[inflow_tag] = {'flux': Constant(-0.22)}
        boundary_conditions[outflow_tag] = {'elev': Constant(0.397)}
        return boundary_conditions

    def set_boundary_conditions_tracer(self, sed_model):
        boundary_conditions = {}
        inflow_tag = 1
        boundary_conditions[inflow_tag] = {'value': sed_model.equiltracer}        
        return boundary_conditions

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

    def get_update_forcings(self, solver_obj):
        return None

    def initialise_fields(self, mesh2d, inputdir, outputdir):
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
    
        DG_2d = get_functionspace(mesh2d, "DG", 1)
        # elevation
        with timed_stage('initialising elevation'):
            chk = DumbCheckpoint(inputdir + "/elevation", mode=FILE_READ)
            elev_init = Function(DG_2d, name="elevation")
            chk.load(elev_init)
            File(outputdir + "/elevation_imported.pvd").write(elev_init)
            chk.close()
        # velocity
        with timed_stage('initialising velocity'):
            chk = DumbCheckpoint(inputdir + "/velocity", mode=FILE_READ)
            V = VectorFunctionSpace(mesh2d, "DG", 1)
            uv_init = Function(V, name="velocity")
            chk.load(uv_init)
            File(outputdir + "/velocity_imported.pvd").write(uv_init)
            chk.close()
        return elev_init, uv_init,


    def get_export_func(self, solver_obj):
        self.bath_export = solver_obj.fields.bathymetry_2d

        def export_func():
            self.bath_file.write(self.bath_export)
        return export_func

    def set_boundary_surface(self):
        """Set the initial displacement of the boundary elevation."""
        pass
