from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time

from adapt_utils.test_cases.trench_test.options import TrenchOptions
from adapt_utils.swe.tsunami.solver import TsunamiProblem

t1 = time.time()

nx = 0.4

op = TrenchOptions(approach='monge_ampere',
                    plot_timeseries=False,
                    plot_pvd=True,
                    debug=False,
                    nonlinear_method='relaxation',
                    # nonlinear_method='quasi_newton',  # FIXME
                    num_adapt=1,
                    qoi_mode='inundation_volume',
                    friction = 'nikuradse',
                    nx=nx,
                    ny = 1,
                    r_adapt_rtol=1.0e-3)

tp = TsunamiProblem(op, levels=0)
tp.setup_solver()


def gradient_interface_monitor(mesh, alpha = 1000.0, beta = 10.0):
    """
    Monitor function focused around the steep_gradient (budd acta numerica)

    NOTE: Defined on the *computational* mesh.

    :kwarg alpha: controls the size of the dense region surrounding the coast.
    """
    P1 = FunctionSpace(mesh, "CG", 1)

    eta = tp.solution.split()[1]
    b = tp.solver_obj.fields.bathymetry_2d
    current_mesh = eta.function_space().mesh()
    P1_current = FunctionSpace(current_mesh, "CG", 1)
    bath_dx_sq = interpolate(pow(b.dx(0), 2), P1_current)
    bath_dy_sq = interpolate(pow(b.dx(1), 2), P1_current)
    bath_dx_dx_sq = interpolate(pow(bath_dx_sq.dx(0), 2), P1_current)
    bath_dy_dy_sq = interpolate(pow(bath_dy_sq.dx(1), 2), P1_current)
    #norm = interpolate(conditional(bath_dx_dx_sq + bath_dy_dy_sq > 10**(-7), bath_dx_dx_sq + bath_dy_dy_sq, Constant(10**(-7))), P1_current)
    norm_two = interpolate(bath_dx_dx_sq + bath_dy_dy_sq, P1_current)
    norm_one = interpolate(bath_dx_sq + bath_dy_sq, P1_current)
    #norm_tmp = interpolate(bath_dx_sq/norm, P1_current)
    norm_one_proj = project(norm_one, P1)
    norm_two_proj = project(norm_two, P1)


    return sqrt(1.0 + (alpha*norm_two_proj) + (beta*norm_one_proj))

tp.monitor_function = gradient_interface_monitor
tp.solve(uses_adjoint=False)

t2 = time.time()



data = pd.read_csv('experimental_data.csv', header = None)

datathetis = []
bathymetrythetis1 = []
diff_thetis = []
for i in range(len(data[0].dropna()[0:3])):
    print(i)
    datathetis.append(data[0].dropna()[i])
    bathymetrythetis1.append(-tp.solver_obj.fields.bathymetry_2d.at([np.round(data[0].dropna()[i],3), 0.55]))
    diff_thetis.append((data[1].dropna()[i] - bathymetrythetis1[-1])**2)
for i in range(4, len(data[0].dropna())):
    print(i)
    datathetis.append(data[0].dropna()[i])
    bathymetrythetis1.append(-tp.solver_obj.fields.bathymetry_2d.at([np.round(data[0].dropna()[i],3), 0.55]))
    diff_thetis.append((data[1].dropna()[i] - bathymetrythetis1[-1])**2)
    
    
df = pd.concat([pd.DataFrame(datathetis, columns = ['x']), pd.DataFrame(bathymetrythetis1, columns = ['bath'])], axis = 1)

df.to_csv('bed_trench_output' + str(nx) + '.csv')

plt.plot(datathetis, bathymetrythetis1, '.', linewidth = 2, label = 'adapted mesh')
plt.legend()
plt.show()
    
print("L2 norm: ")
print(np.sqrt(sum(diff_thetis)))  

print("total time: ")  
print(t2 -t1)

f = open("output_" + str(nx) +'_' + str(1000)+ '.txt', "w+")
f.write(str(np.sqrt(sum(diff_thetis))))
f.write("/n")
f.write(str(t2-t1))
f.close()
