#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:13:47 2019

@author: mc4117
"""

import thetis as th
from thetis import *
import numpy as np
import pandas as pd
import pylab as plt
import time

from adapt_utils.test_cases.trench_test.options import TrenchOptions

timestep = 0.3

t1 = time.time()

op = TrenchOptions(friction = 'nikuradse', nx = 0.2)

op.quadratic_drag_coefficient = Function(op.P1DG).interpolate(op.get_cfactor())

solver_obj = solver2d.FlowSolver2d(op.default_mesh, op.bathymetry)

# Initial conditions

options = solver_obj.options


options.check_volume_conservation_2d = True


# Timestepping
options.timestep = op.dt
options.simulation_export_time = op.dt*op.dt_per_export
options.simulation_end_time = op.end_time
options.timestepper_type = op.timestepper

if hasattr(options.timestepper_options, 'implicitness_theta'):
    options.timestepper_options.implicitness_theta = op.implicitness_theta


if op.solve_tracer:
            options.fields_to_export = ['uv_2d', 'elev_2d', 'tracer_2d'] if op.plot_pvd else []
            options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d', 'tracer_2d'] if op.save_hdf5 else []

# Parameters
options.horizontal_viscosity = Constant(op.base_viscosity)
options.horizontal_diffusivity = Constant(op.base_diffusivity)
options.quadratic_drag_coefficient = op.quadratic_drag_coefficient
options.use_lax_friedrichs_tracer = False 
options.solve_tracer = op.solve_tracer
options.use_tracer_conservative_form = False
options.tracer_advective_velocity_factor = op.corrective_velocity_factor
options.tracer_source_2d = op.source
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1.0

# Boundary conditions
solver_obj.bnd_functions['shallow_water'] = op.set_boundary_conditions(op.P1)
solver_obj.bnd_functions['tracer'] = op.set_boundary_conditions_tracer(op.P1)

tracer_init = op.set_tracer_init(op.P1DG)
solver_obj.assign_initial_conditions(uv=op.uv_init, elev=op.elev_init, tracer=tracer_init)

op.suspended = True
op.bedload = True
op.convective_vel_flag = True

update_forcings_tracer = op.get_update_forcings(solver_obj)

# run model
solver_obj.iterate(update_forcings=update_forcings_tracer)

t2 = time.time()

new_mesh = th.RectangleMesh(16*5*5, 5*1, 16, 1.1)

bath = th.Function(th.FunctionSpace(new_mesh, "CG", 1)).project(solver_obj.fields.bathymetry_2d)

data = pd.read_csv('experimental_data.csv', header=None)
plt.scatter(data[0], data[1], label='Experimental Data')

datathetis = []
bathymetrythetis1 = []
diff_thetis = []
for i in range(len(data[0].dropna())):
    print(i)
    datathetis.append(data[0].dropna()[i])
    bathymetrythetis1.append(-bath.at([np.round(data[0].dropna()[i], 3), 0.55]))
    diff_thetis.append((data[1].dropna()[i] - bathymetrythetis1[-1])**2)

df = pd.concat([pd.DataFrame(datathetis, columns=['x']), pd.DataFrame(bathymetrythetis1, columns=['bath'])], axis=1)

df.to_csv('fixed_output/' + str(nx) + '_bed_trench_output.csv')

print("Time: ")
print(t2-t1)

print("L2 norm: ")
print(np.sqrt(sum(diff_thetis)))

f = open("fixed_output/output_nx" + str(nx) + '_' + str(timestep) + '.txt', "w+")
f.write(str(np.sqrt(sum(diff_thetis))))
f.write("\n")
f.write(str(t2-t1))
f.close()
