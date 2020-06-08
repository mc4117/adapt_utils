from thetis import *

from adapt_utils.test_cases.beach.options import BeachOptions
from adapt_utils.swe.solver import UnsteadyShallowWaterProblem

import time
import datetime

nx = 0.25

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
                   r_adapt_rtol=1.0e-3,
                   init = True,
                   output_dir = outputdir)

swp = UnsteadyShallowWaterProblem(op, levels=0)
swp.setup_solver()


def tracer_interface_monitor(mesh, alpha=1.0, beta=3000.0):
    """
    Monitor function focused around the wet-dry interface.

    NOTE: Defined on the *computational* mesh.

    :kwarg alpha: controls the size of the dense region surrounding the coast.
    :kwarg beta: controls the level of refinement in this region.
    """
    P1 = FunctionSpace(mesh, "DG", 1)
    b = swp.solver_obj.fields.tracer_2d
    #current_mesh = b.function_space().mesh()
    #P1_current = FunctionSpace(current_mesh, "CG", 1)
    diff_proj = project(b, P1)
    
    mon_init = project(1.0 + alpha*pow(cosh(beta*diff_proj), 2), P1)
    
    H = Function(P1)

    tau = TestFunction(P1)
    n = FacetNormal(mesh)
    
    K = 10**6
    a = (inner(tau, H)*dx)+(K*inner(tau.dx(1), H.dx(1))*dx) - inner(tau, mon_init)*dx
    solve(a == 0, H)    

    return H


swp.monitor_function = tracer_interface_monitor
swp.solve(uses_adjoint=False)
