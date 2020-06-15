from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time
import datetime

from adapt_utils.test_cases.beach.options import BeachOptions
from adapt_utils.swe.solver import UnsteadyShallowWaterProblem

t1 = time.time()

nx = 1.5

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
outputdir = 'outputs' + st

op = BeachOptions(approach='fixed_mesh',
                   plot_timeseries=False,
                   plot_pvd=True,
                   debug=False,
                   nonlinear_method='relaxation',
                   num_adapt=1,
                   friction='manning',
                   nx=nx,
                   ny=1,
                   input_dir = 'hydrodynamics_beach_l_sep_nx_330.0',
                   output_dir = outputdir,
                   r_adapt_rtol=1.0e-3,
                   init = True)


swp = UnsteadyShallowWaterProblem(op, levels=0)
swp.setup_solver()

t1 = time.time()

swp.solve(uses_adjoint=False)

t2 = time.time()

print(t2 - t1)

xaxisthetis1 = []
baththetis1 = []

for i in np.linspace(0, 219, 220):
    xaxisthetis1.append(i)
    baththetis1.append(-swp.solver_obj.fields.bathymetry_2d.at([i, 5]))
df = pd.concat([pd.DataFrame(xaxisthetis1, columns = ['x']), pd.DataFrame(baththetis1, columns = ['bath'])], axis = 1)
df.to_csv("final_result_nx" + str(nx) + ".csv", index = False)
