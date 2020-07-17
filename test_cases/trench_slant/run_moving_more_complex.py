
from thetis import *

import firedrake as fire
import pylab as plt
import pandas as pd
import numpy as np
import time
import datetime
from firedrake.petsc import PETSc

from adapt_utils.test_cases.trench_slant.options import TrenchSlantOptions
from adapt_utils.swe.morphological.solver import UnsteadyShallowWaterProblem
from adapt_utils.adapt import recovery
from adapt_utils.norms import local_frobenius_norm, local_norm

t1 = time.time()

nx = 0.2
alpha = 350
beta = 0.5
gamma = 1

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
outputdir = 'outputs' + st

dir = 'hydrodynamics_trench_slant_' + str(nx)
#dir = 'hydrodynamics_trench_slant_0.250.4'
print(dir)

op = TrenchSlantOptions(approach='monge_ampere',
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


def gradient_interface_monitor(mesh, alpha=alpha, beta = beta, gamma=gamma):

    """
    Monitor function focused around the steep_gradient (budd acta numerica)

    NOTE: Defined on the *computational* mesh.

    """
    P1 = FunctionSpace(mesh, "CG", 1)

    b = swp.solver_obj.fields.bathymetry_2d
    bath_gradient = recovery.construct_gradient(b)
    bath_hess = recovery.construct_hessian(b, op=op)
    frob_bath_hess = Function(b.function_space()).project(local_frobenius_norm(bath_hess))
    frob_bath_norm = Function(b.function_space()).project(frob_bath_hess/max(frob_bath_hess.dat.data[:]))
    current_mesh = b.function_space().mesh()
    l2_bath_grad = Function(b.function_space()).project(local_norm(bath_gradient))
    bath_dx_l2_norm = Function(b.function_space()).interpolate(l2_bath_grad/max(l2_bath_grad.dat.data[:]))
    
    comp = interpolate(conditional(beta*bath_dx_l2_norm > gamma*frob_bath_norm, beta*bath_dx_l2_norm, gamma*frob_bath_norm), b.function_space())  
    comp_new = project(comp, P1)
    comp_new2 = interpolate(conditional(comp_new > Constant(0.0), comp_new, Constant(0.0)), P1)
    mon_init = project(sqrt(Constant(1.0) + alpha * comp_new2), P1)

    #K = 10*(0.2**2)/4
    #a = (inner(tau, H)*dx)+(K*inner(grad(tau), grad(H))*dx) - (K*(tau*inner(grad(H), n)))*ds
    #a -= inner(tau, mon_init)*dx
    #solve(a == 0, H)

    return mon_init


swp.monitor_function = gradient_interface_monitor
swp.solve(uses_adjoint=False)

t2 = time.time()

new_mesh = RectangleMesh(16*5*4, 5*4, 16, 1.1)

bath = Function(FunctionSpace(new_mesh, "CG", 1)).project(swp.solver_obj.fields.bathymetry_2d)

def initialise_fields(mesh2d, inputdir):
    """
    Initialise simulation with results from a previous simulation
    """
    V = FunctionSpace(mesh2d, 'CG', 1)
    # elevation
    with timed_stage('initialising bathymetry'):
        chk = DumbCheckpoint(inputdir + "/bathymetry", mode=FILE_READ)
        bath = Function(V, name="bathymetry")
        chk.load(bath)
        chk.close()
        
    return bath

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

export_final_state("adapt_output/hydrodynamics_trench_slant_bath_c"+str(alpha) + "_" + str(beta) + '_' + str(gamma) + '-' + str(nx), bath)

bath_real = initialise_fields(new_mesh, 'hydrodynamics_trench_slant_bath_new_4.0')

print(nx)
print(alpha)
print('L2')
print(fire.errornorm(bath, bath_real))

print("total time: ")
print(t2-t1)

print(beta)
print(gamma)
