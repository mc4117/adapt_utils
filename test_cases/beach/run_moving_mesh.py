from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time

from adapt_utils.test_cases.beach.options import BeachOptions
from adapt_utils.swe.solver import UnsteadyShallowWaterProblem
from adapt_utils.adapt import recovery
from adapt_utils.norms import local_frobenius_norm

nx = 0.25


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
                   r_adapt_rtol=1.0e-3,
                   init = True)

swp = UnsteadyShallowWaterProblem(op, levels=0)
swp.setup_solver()


def wet_dry_interface_monitor(mesh, alpha=20.0, beta=1.0):  
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
    #diff = interpolate(eta + b, P1_current)
    #diff_proj = project(diff, P1)
    horizontal_velocity = interpolate(uv[0], P1_current)    
    #eta_gradient = recovery.construct_gradient(eta)
    #eta_dx_sq = interpolate(pow(eta_gradi
    ent, 2), P1_current)

    uv_gradient = recovery.construct_gradient(horizontal_velocity)
    uv_dx = interpolate(pow(uv_gradient[0], 2), P1_current)
    uv_dx_new = project(uv_dx, P1)
    mon_init = project(sqrt(1.0 + alpha * uv_dx_new), P1)
    #eta_dx_sq_new = project(eta_dx_sq, P1)

    #mon_init = project(sqrt(1.0 + alpha * eta_dx_sq_new), P1)
    #return mon_init
    #return 1.0 + alpha*pow(cosh(beta*diff_proj), -2)
    
    H = Function(P1)
    tau = TestFunction(P1)
    n = FacetNormal(mesh)
    
    K = 10*(0.4**2)/4
    a = (inner(tau, H)*dx)+(K*inner(grad(tau), grad(H))*dx) - (K*(tau*inner(grad(H), n)))*ds
    a -= inner(tau, mon_init)*dx
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

print(sum([(df['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))