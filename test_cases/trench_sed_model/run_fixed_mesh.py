from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time

from adapt_utils.test_cases.trench_sed_model.options import TrenchOptions
from adapt_utils.swe.morphological.solver import UnsteadyShallowWaterProblem

t1 = time.time()

nx = 0.1

dir = 'hydrodynamics_trench_' + str(nx)

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
outputdir = 'outputs' + st

op = TrenchOptions(approach='fixed_mesh',
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

t1 = time.time()

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

df.to_csv('fixed_output/bed_trench_output_uni_' + str(nx) + '.csv')

plt.plot(datathetis, bathymetrythetis1, '.', linewidth=2, label='fixed mesh')
plt.legend()
plt.show()

#print("L2 norm: ")
#print(np.sqrt(sum(diff_thetis)))
print(nx)
print("total time: ")
print(t2-t1)
