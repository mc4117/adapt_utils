#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 16:37:44 2020

@author: mc4117
"""

import thetis as th
import firedrake as fire
import pylab as plt
import numpy as np

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

bath1 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_0_0.1')
bath2 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_5_0.1')
bath3 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_10_0.1')
bath4 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_15_0.1')
bath5 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_20_0.1')
bath6 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_25_0.1')
bath7 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_30_0.1')
bath8 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_35_0.1')
bath9 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_40_0.1')
bath10 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_45_0.1')
bath11 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_50_0.1')
bath12 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_55_0.1')
bath13 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_60_0.1')
bath14 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_70_0.1')
bath15 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_80_0.1')
bath16 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_90_0.1')
bath17 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_100_0.1')
bath18 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_110_0.1')
bath19 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_120_0.1')
bath20 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_125_0.1')
bath21 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_150_0.1')
bath22 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_175_0.1')
bath23 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_250_0.1')
bath24 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_275_0.1')
bath25 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_300_0.1')
bath26 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_325_0.1')
bath27 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_350_0.1')
bath28 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_400_0.1')
bath_real = initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_new_4.0')

errorlist = []
errorlist.append(fire.errornorm(bath1, bath_real))
errorlist.append(fire.errornorm(bath2, bath_real))
errorlist.append(fire.errornorm(bath3, bath_real))
errorlist.append(fire.errornorm(bath4, bath_real))
errorlist.append(fire.errornorm(bath5, bath_real))
errorlist.append(fire.errornorm(bath6, bath_real))
errorlist.append(fire.errornorm(bath7, bath_real))
errorlist.append(fire.errornorm(bath8, bath_real))
errorlist.append(fire.errornorm(bath9, bath_real))
errorlist.append(fire.errornorm(bath10, bath_real))
errorlist.append(fire.errornorm(bath11, bath_real))
errorlist.append(fire.errornorm(bath12, bath_real))
errorlist.append(fire.errornorm(bath13, bath_real))
errorlist.append(fire.errornorm(bath14, bath_real))
errorlist.append(fire.errornorm(bath15, bath_real))
errorlist.append(fire.errornorm(bath16, bath_real))
errorlist.append(fire.errornorm(bath17, bath_real))
errorlist.append(fire.errornorm(bath18, bath_real))
errorlist.append(fire.errornorm(bath19, bath_real))
errorlist.append(fire.errornorm(bath20, bath_real))
errorlist.append(fire.errornorm(bath21, bath_real))
errorlist.append(fire.errornorm(bath22, bath_real))
errorlist.append(fire.errornorm(bath23, bath_real))
errorlist.append(fire.errornorm(bath24, bath_real))
errorlist.append(fire.errornorm(bath25, bath_real))
errorlist.append(fire.errornorm(bath26, bath_real))
errorlist.append(fire.errornorm(bath27, bath_real))
errorlist.append(fire.errornorm(bath28, bath_real))


print(errorlist)

bath1 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_0_0.2')
bath2 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_5_0.2')
bath3 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_10_0.2')
bath4 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_15_0.2')
bath5 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_20_0.2')
bath6 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_25_0.2')
bath7 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_30_0.2')
bath8 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_35_0.2')
bath9 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_40_0.2')

errorlist = []
errorlist.append(fire.errornorm(bath1, bath_real))
errorlist.append(fire.errornorm(bath2, bath_real))
errorlist.append(fire.errornorm(bath3, bath_real))
errorlist.append(fire.errornorm(bath4, bath_real))
errorlist.append(fire.errornorm(bath5, bath_real))
errorlist.append(fire.errornorm(bath6, bath_real))
errorlist.append(fire.errornorm(bath7, bath_real))
errorlist.append(fire.errornorm(bath8, bath_real))
errorlist.append(fire.errornorm(bath9, bath_real))

print(errorlist)

bath1 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c100_0.1')
bath2 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c150_0.1')
bath3 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c200_0.1')
bath4 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c250_0.1')
bath5 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c300_0.1')
bath6 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c350_0.1')

errorlist = []
errorlist.append(fire.errornorm(bath1, bath_real))
errorlist.append(fire.errornorm(bath2, bath_real))
errorlist.append(fire.errornorm(bath3, bath_real))
errorlist.append(fire.errornorm(bath4, bath_real))
errorlist.append(fire.errornorm(bath5, bath_real))
errorlist.append(fire.errornorm(bath6, bath_real))

print(errorlist)

bath1 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_0_0.4')
bath2 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_1_0.4')
bath3 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_2_0.4')
bath4 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_3_0.4')
bath5 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_4_0.4')
#bath6 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_5_0.4')
bath7 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_7_0.4')
bath8 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_9_0.4')
bath9 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_10_0.4')
bath10 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_15_0.4')
bath11 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_20_0.4')

errorlist = []
errorlist.append(fire.errornorm(bath1, bath_real))
errorlist.append(fire.errornorm(bath2, bath_real))
errorlist.append(fire.errornorm(bath3, bath_real))
errorlist.append(fire.errornorm(bath4, bath_real))
errorlist.append(fire.errornorm(bath5, bath_real))
#errorlist.append(fire.errornorm(bath6, bath_real))
errorlist.append(fire.errornorm(bath7, bath_real))
errorlist.append(fire.errornorm(bath8, bath_real))
errorlist.append(fire.errornorm(bath9, bath_real))
errorlist.append(fire.errornorm(bath10, bath_real))
errorlist.append(fire.errornorm(bath11, bath_real))

print(errorlist)

bath1 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c0_0.4')
bath2 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c5_0.4')
bath3 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c10_0.4')
bath4 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c15_0.4')
bath5 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c20_0.4')
bath6 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c25_0.4')
bath7 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c30_0.4')
bath8 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c35_0.4')
bath9 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c40_0.4')
bath10 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c45_0.4')
bath11 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c50_0.4')
bath12 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c55_0.4')
bath13 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c60_0.4')
bath14 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c65_0.4')
bath15 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c70_0.4')
bath16 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c75_0.4')
bath17 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c80_0.4')
bath18 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c85_0.4')
bath19 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c90_0.4')
bath20 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c100_0.4')
bath21 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c120_0.4')
bath22 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c150_0.4')
bath23 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c175_0.4')
bath24 = initialise_fields(new_mesh, 'adapt_output/hydrodynamics_trench_slant_bath_c200_0.4')

errorlist = []
errorlist.append(fire.errornorm(bath1, bath_real))
errorlist.append(fire.errornorm(bath2, bath_real))
errorlist.append(fire.errornorm(bath3, bath_real))
errorlist.append(fire.errornorm(bath4, bath_real))
errorlist.append(fire.errornorm(bath5, bath_real))
errorlist.append(fire.errornorm(bath6, bath_real))
errorlist.append(fire.errornorm(bath7, bath_real))
errorlist.append(fire.errornorm(bath8, bath_real))
errorlist.append(fire.errornorm(bath9, bath_real))
errorlist.append(fire.errornorm(bath10, bath_real))
errorlist.append(fire.errornorm(bath11, bath_real))
errorlist.append(fire.errornorm(bath12, bath_real))
errorlist.append(fire.errornorm(bath13, bath_real))
errorlist.append(fire.errornorm(bath14, bath_real))
errorlist.append(fire.errornorm(bath15, bath_real))
errorlist.append(fire.errornorm(bath16, bath_real))
errorlist.append(fire.errornorm(bath17, bath_real))
errorlist.append(fire.errornorm(bath18, bath_real))
errorlist.append(fire.errornorm(bath19, bath_real))
errorlist.append(fire.errornorm(bath20, bath_real))
errorlist.append(fire.errornorm(bath21, bath_real))
errorlist.append(fire.errornorm(bath22, bath_real))
errorlist.append(fire.errornorm(bath23, bath_real))
errorlist.append(fire.errornorm(bath24, bath_real))

print(errorlist)