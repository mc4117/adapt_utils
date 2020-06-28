#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 16:37:44 2020

@author: mc4117
"""

import thetis as th
import firedrake as fire

def initialise_fields(mesh2d, inputdir):
    """
    Initialise simulation with results from a previous simulation
    """
    V = th.FunctionSpace(mesh2d, 'CG', 1)
    # elevation
    with th.timed_stage('initialising bathymetry'):
        chk = th.DumbCheckpoint(inputdir + "/bathymetry", mode=th.FILE_READ)
        bath = th.Function(V, name="bathymetry")
        chk.load(bath)
        chk.close()
        
    return bath

new_mesh = th.RectangleMesh(16*5*4, 5*4, 16, 1.1)

bath1 = initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_0.1')
bath2 = initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_0.2')
bath3 = initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_0.4')
bath4 = initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_0.8')
bath5 =  initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_2.0')
bath_real = initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_4.0')

errorlist = []
errorlist.append(fire.errornorm(bath1, bath_real))
errorlist.append(fire.errornorm(bath2, bath_real))
errorlist.append(fire.errornorm(bath3, bath_real))
errorlist.append(fire.errornorm(bath4, bath_real))
errorlist.append(fire.errornorm(bath5, bath_real))

print(errorlist)
