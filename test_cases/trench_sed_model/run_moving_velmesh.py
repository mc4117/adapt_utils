from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time
import datetime

from adapt_utils.test_cases.trench_sed_model.options import TrenchOptions
from adapt_utils.swe.morphological.solver import UnsteadyShallowWaterProblem
from adapt_utils.adapt import recovery
from adapt_utils.norms import local_frobenius_norm

t1 = time.time()

nx = 0.4
alpha = 0.15

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
outputdir = 'outputs' + st

dir = 'hydrodynamics_trench_' + str(nx)

op = TrenchOptions(approach='monge_ampere',
                   input_dir = dir,
                   output_dir = outputdir,                   
                   plot_timeseries=False,
                   plot_pvd=True,
                   debug=False,
                   nonlinear_method='relaxation',
                   num_adapt=1,
                   friction='nikuradse',
                   nx=nx,
                   ny=1,
                   r_adapt_rtol=1.0e-3,
                   init = True)

swp = UnsteadyShallowWaterProblem(op, levels=0)
swp.setup_solver()


def vel_interface_monitor(mesh, alpha=alpha, beta=1.0):
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
    
    abs_horizontal_velocity = interpolate(abs(horizontal_velocity), P1_current)


    uv_gradient = recovery.construct_gradient(horizontal_velocity)
    div_uv = interpolate(sqrt(inner(uv_gradient, uv_gradient)), P1_current)
    div_uv_star = interpolate(conditional(div_uv/(beta*max(div_uv.dat.data[:])) < Constant(1), 
                                          div_uv/(beta*max(div_uv.dat.data[:])) , Constant(1)), P1_current)
    
    abs_uv_star = interpolate(conditional(abs_horizontal_velocity/(beta*max(abs_horizontal_velocity.dat.data[:])) < Constant(1), 
                                     abs_horizontal_velocity/(beta*max(abs_horizontal_velocity.dat.data[:])) , Constant(1)), P1_current)
    
    comp = interpolate(conditional(abs_uv_star > div_uv_star, abs_uv_star, div_uv_star)**2, P1_current)      
    comp_new = project(comp, P1)
    comp_new2 = interpolate(conditional(comp_new > Constant(0.0), comp_new, Constant(0.0)), P1)
    mon_init = project(sqrt(1.0 + alpha * comp_new2), P1)
 

    return mon_init


swp.monitor_function = vel_interface_monitor
swp.solve(uses_adjoint=False)

t2 = time.time()

new_mesh = RectangleMesh(16*5*5, 5*1, 16, 1.1)

bath = Function(FunctionSpace(new_mesh, "CG", 1)).project(swp.solver_obj.fields.bathymetry_2d)

data = pd.read_csv('experimental_data.csv', header=None)

datathetis = []
bathymetrythetis1 = []
diff_thetis = []
for i in np.linspace(0, 15.9, 160):
    datathetis.append(i)
    bathymetrythetis1.append(-bath.at([i, 0.55]))
    #diff_thetis.append((data[1].dropna()[i] - bathymetrythetis1[-1])**2)

df = pd.concat([pd.DataFrame(datathetis, columns=['x']), pd.DataFrame(bathymetrythetis1, columns=['bath'])], axis=1)

#df.to_csv('adapt_output/bed_trench_output_uni_' + str(nx) + '_' + str(alpha) + '.csv')

datathetis = []
bathymetrythetis1 = []
diff_thetis = []
for i in range(len(data[0].dropna())):
    print(i)
    datathetis.append(data[0].dropna()[i])
    bathymetrythetis1.append(-bath.at([np.round(data[0].dropna()[i], 3), 0.55]))
    diff_thetis.append((data[1].dropna()[i] - bathymetrythetis1[-1])**2)

df_exp = pd.concat([pd.DataFrame(datathetis, columns=['x']), pd.DataFrame(bathymetrythetis1, columns=['bath'])], axis=1)

#df_exp.to_csv('adapt_output/bed_trench_output_' + str(nx) + '_' + str(alpha) + '.csv')



print("Total error: ")
print(np.sqrt(sum(diff_thetis)))

print("total time: ")
print(t2-t1)


df_real = pd.read_csv('fixed_output/bed_trench_output_uni_4.csv')
print("Mesh error: ")
print(sum([(df['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
#f = open("adapt_output/output_frob_norm_" + str(nx) + '_' + str(alpha) + '.txt', "w+")
#f.write(str(np.sqrt(sum(diff_thetis))))
#f.write("\n")
#f.write(str(t2-t1))
#f.close()
