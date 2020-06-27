from thetis import *

import pylab as plt
import pandas as pd
import numpy as np
import time
from firedrake.petsc import PETSc

from adapt_utils.test_cases.trench_slant.options import TrenchSlantOptions
from adapt_utils.swe.morphological.solver import UnsteadyShallowWaterProblem

def export_final_state(inputdir, bathymetry_2d):
    """
    Export fields to be used in a subsequent simulation
    """
    if not os.path.exists(inputdir):
        os.makedirs(inputdir)
    print_output("Exporting fields for subsequent simulation")

    chk = DumbCheckpoint(inputdir + "/bathymetry", mode=FILE_CREATE)
    chk.store(bathymetry_2d, name="bathymetry")
    File(inputdir + '/bathout.pvd').write(bathymetry_2d)
    chk.close()
    
    plex = bathymetry_2d.function_space().mesh()._plex
    viewer = PETSc.Viewer().createHDF5(inputdir + '/myplex.h5', 'w')
    viewer(plex)        

t1 = time.time()

nx = 4.0

dir = 'hydrodynamics_trench_slant_' + str(nx)

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
outputdir = 'outputs' + st

op = TrenchSlantOptions(approach='fixed_mesh',
                   input_dir = dir,
                   output_dir = outputdir,
                   plot_timeseries=False,
                   plot_pvd=True,
                   debug=False,
                   nonlinear_method='relaxation',
                   num_adapt=1,
                   friction='nikuradse',
                   nx=nx,
                   ny=nx,
                   r_adapt_rtol=1.0e-3,
                   init = True)

swp = UnsteadyShallowWaterProblem(op, levels=0)
swp.setup_solver()

t1 = time.time()

swp.solve(uses_adjoint=False)

t2 = time.time()

new_mesh = RectangleMesh(16*5*4, 5*4, 16, 1.1)

bath = Function(FunctionSpace(new_mesh, "CG", 1)).project(swp.solver_obj.fields.bathymetry_2d)

export_final_state("hydrodynamics_trench_slant_bath_"+str(nx), bath)

print("total time: ")
print(t2-t1)

