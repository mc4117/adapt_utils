from thetis import *

from adapt_utils.test_cases.inundated_beach.options import BalzanoOptions
from adapt_utils.swe.tsunami.solver import TsunamiProblem

op = BalzanoOptions(approach='monge_ampere',
                    plot_timeseries=False,  # FIXME
                    plot_pvd=True,
                    debug=True,
                    nonlinear_method='relaxation',
                    # nonlinear_method='quasi_newton',  # FIXME
                    num_adapt=1,
                    qoi_mode='inundation_volume',
                    n=2,
                    r_adapt_rtol=1.0e-3)
<<<<<<< HEAD

#op.default_mesh = RectangleMesh(17, 2, 1.5*op.basin_x, 1200.0)

swp = UnsteadyShallowWaterProblem(op, levels=0)
swp.setup_solver()
=======
tp = TsunamiProblem(op, levels=0)
tp.setup_solver()
>>>>>>> upstream/master


def wet_dry_interface_monitor(mesh, alpha=1.0, beta=1.0):  # FIXME: all this projection is expensive!
    """
    Monitor function focused around the wet-dry interface.

    NOTE: Defined on the *computational* mesh.

    :kwarg alpha: controls the size of the dense region surrounding the coast.
    :kwarg beta: controls the level of refinement in this region.
    """
    P1 = FunctionSpace(mesh, "CG", 1)
    eta = swp.solution.split()[1]
    b = swp.fields['bathymetry']
    current_mesh = eta.function_space().mesh()
    P1_current = FunctionSpace(current_mesh, "CG", 1)
    diff = interpolate(eta + b, P1_current)
    diff_proj = project(diff, P1)
    return 1.0 + alpha*pow(cosh(beta*diff_proj), -2)

tp.monitor_function = wet_dry_interface_monitor
tp.solve(uses_adjoint=False)
# TODO: Evaluate QoI properly, accounting for mesh adaptation
