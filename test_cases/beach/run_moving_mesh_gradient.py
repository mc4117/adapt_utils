from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time
import datetime

from adapt_utils.test_cases.beach.options import BeachOptions
from adapt_utils.swe.solver import UnsteadyShallowWaterProblem
from adapt_utils.adapt import recovery
from adapt_utils.norms import local_frobenius_norm

nx = 0.25

alpha_star = 20.0

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
outputdir = 'outputs' + st

op = BeachOptions(approach='monge_ampere',
                   plot_timeseries=False,
                   plot_pvd=True,
                   debug=False,
                   nonlinear_method='relaxation',
                   num_adapt=1,
                   friction='manning',
                   nx=nx,
                   ny=1,
                   input_dir = 'hydrodynamics_beach_l_sep_nx_55.0',
                   output_dir = outputdir,
                   r_adapt_rtol=1.0e-3,
                   init = True)

swp = UnsteadyShallowWaterProblem(op, levels=0)
swp.setup_solver()


def wet_dry_interface_monitor(mesh, alpha=alpha_star, beta=1.0):
    """
    Monitor function focused around the wet-dry interface.

    NOTE: Defined on the *computational* mesh.

    :kwarg alpha: controls the size of the dense region surrounding the coast.
    :kwarg beta: controls the level of refinement in this region.
    """
    P1 = FunctionSpace(mesh, "CG", 1)
    eta = swp.solution.split()[1]
    uv = swp.solution.split()[0]


    #b = swp.solver_obj.fields.bathymetry_2d
    current_mesh = eta.function_space().mesh()
    P1_current = FunctionSpace(current_mesh, "CG", 1)

    horizontal_velocity = interpolate(uv[0], P1_current)    


    uv_gradient = recovery.construct_gradient(horizontal_velocity)
    uv_dx = interpolate(pow(uv_gradient[0], 2), P1_current)
    uv_dy = interpolate(pow(uv_gradient[1], 2), P1_current)
    div_uv = interpolate(sqrt(uv_dx + uv_dy), P1_current)
    div_uv_star = interpolate(conditional(div_uv/(beta*max(div_uv.dat.data[:])) < Constant(1), 
                                          div_uv/(beta*max(div_uv.dat.data[:])) , Constant(1)), P1_current)

    comp = interpolate(div_uv_star**2, P1_current)
    comp_new = project(comp, P1)
    comp_new2 = interpolate(conditional(comp_new > Constant(0.0), comp_new, Constant(0.0)), P1)
    mon_init = project(sqrt(1.0 + alpha * comp_new2), P1)

    H = Function(P1)
    tau = TestFunction(P1)
    
    K = 100*(0.4**2)/4
    a = (inner(tau, H)*dx)+(K*inner(tau.dx(1), H.dx(1))*dx) - inner(tau, mon_init)*dx
    solve(a == 0, H)    

    return H


swp.monitor_function = wet_dry_interface_monitor
swp.solve(uses_adjoint=False)

t2 = time.time()

new_mesh = RectangleMesh(220*2, 10, 220, 10)

bath = Function(FunctionSpace(new_mesh, "CG", 1)).project(swp.solver_obj.fields.bathymetry_2d)

xaxisthetis1 = []
baththetis1 = []

for i in np.linspace(0, 219, 220):
    xaxisthetis1.append(i)
    baththetis1.append(-bath.at([i, 5]))
    
df = pd.concat([pd.DataFrame(xaxisthetis1, columns = ['x']), pd.DataFrame(baththetis1, columns = ['bath'])], axis = 1)

df_real = pd.read_csv('final_result_nx2.csv')

print(alpha_star)
print(sum([(df['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))

print("total time: ")
print(t2-t1)

f = open("adapt_output/output_grad_norm_" + str(nx) + '_' + str(alpha_star) + '.txt', "w+")
f.write(str(sum([(df['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))])))
f.write("\n")
f.write(str(t2-t1))
f.close()