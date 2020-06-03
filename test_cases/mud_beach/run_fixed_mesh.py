from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time

from adapt_utils.test_cases.mud_beach.options import MudBeachOptions
from adapt_utils.swe.solver import UnsteadyShallowWaterProblem

t1 = time.time()

nx = 0.25

op = MudBeachOptions(approach='fixed_mesh',
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

print(t2 - t1)

xaxisthetis1 = []
baththetis1 = []

for i in np.linspace(0, 9900, 9901):
    xaxisthetis1.append(i)
    baththetis1.append(-swp.solver_obj.fields.bathymetry_2d.at([i, 400]))
df = pd.concat([pd.DataFrame(xaxisthetis1, columns = ['x']), pd.DataFrame(baththetis1, columns = ['bath'])], axis = 1)
df.to_csv("final_result_nx_2_" + str(nx) + ".csv", index = False)
