#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:13:47 2019

@author: mc4117
"""

import thetis as th
import morphological_hydro_mud_source_sink as morph
import numpy as np
import pandas as pd
import pylab as plt
import os

def boundary_conditions_fn_balzano(bathymetry_2d, flag = None, morfac = 1, t_new = 0, state = 'initial'):
    """
    Define boundary conditions for problem to be used in morphological section.
    
    Inputs:
    morfac - morphological scale factor used when calculating time dependent boundary conditions
    t_new - timestep model currently at used when calculating time dependent boundary conditions
    state - when 'initial' this is the initial boundary condition set; when 'update' these are the boundary
            conditions set during update forcings (ie. if fluc_bcs = True, this will be called)
    """
    left_bnd_id = 1
    right_bnd_id = 3
    left_string = ['flux']#['elev']
    right_string = []

    # set boundary conditions
    swe_bnd = {}
    
    # boundary conditions
    h_amp = -2      # ocean boundary forcing amplitude
    v_amp = -0.4
    h_T = 12*3600.    # ocean boundary forcing period
    ocean_elev_func = lambda t: h_amp * th.sin(2 * th.pi * t / h_T)
    ocean_vel_func = lambda t: v_amp * th.cos(2 * th.pi * t / h_T)
    
    if state == 'initial':

        ocean_elev = ocean_elev_func(0.0)
        ocean_vel = ocean_vel_func(0.0)
        ocean_flux = -ocean_vel*(ocean_elev+4.5)*ly
        
        print(ocean_flux)

        inflow_constant = [ocean_flux]#[ocean_elev]
        outflow_constant = []
        return swe_bnd, left_bnd_id, right_bnd_id, inflow_constant, outflow_constant, left_string, right_string
    elif state == 'update':
        ocean_elev = ocean_elev_func(t_new)
        ocean_vel = ocean_vel_func(t_new)
        depth = th.Function(P1_2d).interpolate(solver_obj.depth.get_total_depth(ocean_elev))
        ocean_flux = -ocean_vel*(depth.at([0.0, 400.0]))*ly
        print(ocean_flux)
        inflow_constant = [ocean_flux]
        outflow_constant = []

        return inflow_constant, outflow_constant
# define mesh
#mesh2d = th.Mesh('extr_mesh.msh')
lx = 10000
ly = 800
nx = 50
ny = 4
#mesh2d = th.Mesh('extr_mesh_2.msh')
mesh2d = th.RectangleMesh(nx, ny, lx, ly)

x,y = th.SpatialCoordinate(mesh2d)

# define function spaces
V = th.FunctionSpace(mesh2d, 'CG', 1)
P1_2d = th.FunctionSpace(mesh2d, 'DG', 1)
vectorP1_2d = th.VectorFunctionSpace(mesh2d, 'DG', 1)


# define underlying bathymetry
bathymetry_2d = th.Function(V, name='Bathymetry')


bathymetry_2d.interpolate(4.5 - x/1000)

# define initial elevation
elev_init = th.Function(P1_2d).interpolate(th.Constant(0.0))
uv_init = th.as_vector((10**(-7), 0.0))

wd_fn = th.Constant(0.2)

solver_obj, update_forcings_hydrodynamics, outputdir = morph.hydrodynamics_only(boundary_conditions_fn_balzano, mesh2d, bathymetry_2d, uv_init, elev_init, wetting_and_drying = True, wetting_alpha = wd_fn, fluc_bcs = True, average_size = 250 * (10**(-6)), dt=125, t_end=24*3600.*10)

# User-defined output: moving bathymetry and eta_tilde
wd_bathfile = th.File(os.path.join(outputdir, 'moving_bath.pvd'))
moving_bath = th.Function(P1_2d, name="moving_bath")
eta_tildefile = th.File(os.path.join(outputdir, 'eta_tilde.pvd'))
eta_tilde = th.Function(P1_2d, name="eta_tilde")

# user-specified export function
def export_func():
    wd_bath_displacement = solver_obj.eq_sw.bathymetry_displacement_mass_term.wd_bathymetry_displacement
    eta = solver_obj.fields.elev_2d
    moving_bath.project(bathymetry_2d + wd_bath_displacement(eta))
    wd_bathfile.write(moving_bath)
    eta_tilde.project(eta+wd_bath_displacement(eta))
    eta_tildefile.write(eta_tilde)

"""
# run hydro model
solver_obj.iterate(update_forcings = update_forcings_hydrodynamics, export_func=export_func)
"""

uv, elev = solver_obj.fields.solution_2d.split()
morph.export_final_state("mud_init", uv, elev)


solver_obj, update_forcings_tracer, diff_bathy, diff_bathy_file, outputdir = morph.morphological(boundary_conditions_fn = boundary_conditions_fn_balzano, morfac = 36.5, morfac_transport = True, suspendedload = True, convectivevel = True,\
                    bedload = True, angle_correction = True, slope_eff = True, seccurrent = False, sediment_slide = False, fluc_bcs = True, diffusivity =0.15, wetting_and_drying = True, wetting_alpha = wd_fn, viscosity_hydro = 10**(-6),\
                    mesh2d = mesh2d, bathymetry_2d = bathymetry_2d, input_dir = 'mud_init', ks = 0.025, average_size = 150*(10**(-6)), dt = 150, final_time = 3650*24*3600, update_bedlevel = True, cons_tracer = True, tracer_init = th.Constant(0.0),  depth_integrated = True)


# run model
solver_obj.iterate(update_forcings = update_forcings_tracer)