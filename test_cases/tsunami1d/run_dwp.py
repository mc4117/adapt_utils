import firedrake

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from adapt_utils.swe.spacetime.solver import SpaceTimeShallowWaterProblem
from adapt_utils.test_cases.tsunami1d.options import Tsunami1dOptions


matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text', usetex=True)

# Parameters
debug = True
plot_pdf = False
plot_pvd = True
save_hdf5 = True
forward = True
adjoint = True

# Spatial discretisation
n = 100
# n = 2000  # (Value used in original paper)
dx = 1/n

# Time discretisation
celerity = 20.0*np.sqrt(9.81)
# dt = 2000.0*dx/celerity
dt = 5.0
# dt = 1.0  # (Value used in original paper)

# NOTE: Forward and adjoint relatively stable with n = 500 and dt = 1.5
op = Tsunami1dOptions(debug=debug, approach='dwp', nx=n, dt=dt,
                      save_hdf5=save_hdf5, plot_pvd=plot_pvd,
                      horizontal_length_scale=1000.0, time_scale=10.0)
op.h_min = 100.0/op.L
op.h_max = 100.0e+3/op.L
op.target = 10000.0
op.num_adapt = 1
# op.norm_order = 1
# op.normalisation = 'error'

swp = SpaceTimeShallowWaterProblem(op, discrete_adjoint=False)
swp.setup_solver_forward()
swp.solve_forward()
swp.setup_solver_adjoint()
swp.solve_adjoint()
swp.dwp_indication()
swp.indicator.interpolate(abs(swp.indicator))
swp.get_isotropic_metric()
swp.adapt_mesh()

v = 0.02
L = op.L  # Horizontal length scale
T = op.T  # Time scale

# FIXME: Solution of equations on new mesh

if forward:
    # Solve forward problem
    swp.setup_solver_forward()
    swp.solve_forward()
    eta = swp.solution.split()[1]

    # Plot forward
    fig = plt.figure(figsize=(3.2, 4.8))
    ax = fig.add_subplot(111)
    # firedrake.plot(firedrake.interpolate(abs(eta), eta.function_space()), axes=ax, vmin=v, vmax=v+0.0001, cmap=matplotlib.cm.Reds)
    firedrake.plot(firedrake.interpolate(abs(eta), eta.function_space()), axes=ax, cmap=matplotlib.cm.Reds)
    ax.invert_xaxis()
    ymin = op.start_time
    ymax = op.end_time
    plt.xlabel("Kilometres offshore")
    plt.ylabel("Hours")
    plt.tight_layout()
    plt.axvline(50e+3/L, ymin=ymin, ymax=ymax, linestyle='--', color='k')
    plt.xlim([400e+3/L, 0.0/L])
    plt.ylim([ymin, ymax])
    plt.xticks([50e+3/L, 150e+3/L, 250e+3/L, 350e+3/L], ["50", "150", "250", "350"])
    plt.yticks([1800.0/T, 3600.0/T], ["0.5", "1.0"])
    fname = os.path.join(op.di, "forward_{:d}".format(n))
    plt.savefig(fname + ".png")
    if plot_pdf:
        plt.savefig(fname + ".pdf")

if adjoint:
    # Solve adjoint problem
    swp.setup_solver_adjoint()
    swp.solve_adjoint()
    zeta = swp.adjoint_solution.split()[1]

    # Plot adjoint
    fig = plt.figure(figsize=(3.2, 4.8))
    ax = fig.add_subplot(111)
    # firedrake.plot(firedrake.interpolate(abs(zeta), zeta.function_space()), axes=ax, vmin=v, vmax=v+0.0001, cmap=matplotlib.cm.Blues)
    firedrake.plot(firedrake.interpolate(abs(zeta), zeta.function_space()), axes=ax, cmap=matplotlib.cm.Blues)
    ax.invert_xaxis()
    ymin = op.start_time
    ymax = op.end_time
    plt.xlabel("Kilometres offshore")
    plt.ylabel("Hours")
    plt.tight_layout()
    plt.axvline(50e+3/L, ymin=ymin, ymax=ymax, linestyle='--', color='k')
    plt.xlim([400e+3/L, 0.0/L])
    plt.ylim([ymin, ymax])
    plt.xticks([50e+3/L, 150e+3/L, 250e+3/L, 350e+3/L], ["50", "150", "250", "350"])
    plt.yticks([], [])
    fname = os.path.join(op.di, "adjoint_{:d}".format(n))
    plt.savefig(fname + ".png")
    if plot_pdf:
        plt.savefig(fname + ".pdf")

if forward and adjoint:
    # Take inner product of forward and adjoint solutions
    swp.dwp_indication()
    swp.plot()
    dwp = swp.indicator

    # Plot inner product
    fig = plt.figure(figsize=(3.2, 4.8))
    ax = fig.add_subplot(111)
    # firedrake.plot(firedrake.interpolate(abs(dwp), dwp.function_space()), axes=ax, vmin=v, vmax=v+0.0001, cmap=matplotlib.cm.Greens)
    firedrake.plot(firedrake.interpolate(abs(dwp), dwp.function_space()), axes=ax, cmap=matplotlib.cm.Greens)
    ax.invert_xaxis()
    ymin = op.start_time
    ymax = op.end_time
    plt.xlabel("Kilometres offshore")
    plt.ylabel("Hours")
    plt.tight_layout()
    plt.axvline(50e+3/L, ymin=ymin, ymax=ymax, linestyle='--', color='k')
    plt.xlim([400e+3/L, 0.0/L])
    plt.ylim([ymin, ymax])
    plt.xticks([50e+3/L, 150e+3/L, 250e+3/L, 350e+3/L], ["50", "150", "250", "350"])
    plt.yticks([], [])
    fname = os.path.join(op.di, "dwp_{:d}".format(n))
    plt.savefig(fname + ".png")
    if plot_pdf:
        plt.savefig(fname + ".pdf")
