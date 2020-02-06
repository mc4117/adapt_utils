from thetis import *

from adapt_utils.case_studies.tohoku.options import TohokuOptions
from adapt_utils.swe.tsunami.solver import TsunamiProblem

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-num_initial_adapt")
parser.add_argument("-n")
parser.add_argument("-end_time")
parser.add_argument("-adapt_field")
parser.add_argument("-approach")
parser.add_argument("-target")
args = parser.parse_args()

n = int(args.n or 40)
num_adapt = int(args.num_initial_adapt or 0)
offset = 0  # TODO: Use offset

op = TohokuOptions(utm=True, plot_pvd=True, n=n, offset=offset)
op.end_time = float(args.end_time or op.end_time)

# Set wetting and drying parameter
# op.wetting_and_drying_alpha.assign(0.5)
h = CellSize(op.default_mesh)
b = op.bathymetry
P0 = FunctionSpace(op.default_mesh, "DG", 0)  # NOTE: alpha is enormous in this approach (O(km))
op.wetting_and_drying_alpha = interpolate(h*sqrt(dot(grad(b), grad(b))), P0)

# Adapt mesh in lonlat space
ext = None
if num_adapt > 0:
    op_init = TohokuOptions(utm=False, n=n, offset=offset, nonlinear_method='relaxation')
    op_init.num_adapt = num_adapt
    op_init.adapt_field = args.adapt_field or 'bathymetry'
    op_init.h_max = 1.0e+10
    op_init.target = float(args.target or 1.0e+4)
    swp_init = TsunamiProblem(op_init, levels=0)
    swp_init.initialise_mesh(approach=args.approach or 'hessian', alpha=0.001, beta=0.005)
    ext = "{:s}_{:d}".format(op_init.adapt_field, op_init.num_adapt)

    # Establish coordinate transformation between lonlat and utm meshes
    op_init.get_utm_mesh()
    op.__init__(mesh=op_init.default_mesh, utm=True, plot_pvd=True, n=n, offset=offset)
    op.get_lonlat_mesh()

swp = TsunamiProblem(op, levels=0)
swp.solve()
op.plot_timeseries("P02", extension=ext)
op.plot_timeseries("P06", extension=ext)
