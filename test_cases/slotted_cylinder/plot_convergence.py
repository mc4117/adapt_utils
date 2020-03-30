import os
import matplotlib.pyplot as plt
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-approach", help="Choose adaptive approach")
parser.add_argument("-mesh_type", help="Choose from 'circle' or 'square'")
args = parser.parse_args()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fontsize = 18

mesh_type = args.mesh_type or 'circle'
approach = args.approach or 'fixed_mesh'

shapes = ('Gaussian', 'Cone', 'Slotted cylinder')

# Load data from log files
elements = []
dat = {}
for shape in shapes:
    dat[shape] = {'qois': [], 'errors': []}
for n in (1, 2, 4, 8):
    with open(os.path.join('outputs', approach, '{:s}_{:d}.log'.format(mesh_type, n)), 'r') as f:
        f.readline()
        f.readline()
        f.readline()
        for shape in shapes:
            line = f.readline().split()
            approx = float(line[-2])
            exact = float(line[-4])
            dat[shape]['qois'].append(approx)
            dat[shape]['errors'].append(100.0*abs(1.0 - approx/exact))
        f.readline()
        elements.append(int(f.readline().split()[-1]))

# for shape in shapes:
#     fig, ax = plt.subplots(figsize=(5, 5))
#     ax.semilogx(elements, dat[shape]['qois'], label=shape)
#     ax.set_xlabel('Element count')
#     ax.set_ylabel('Volume of solid body')

fig, ax = plt.subplots(figsize=(5, 5))
for shape in shapes:
    ax.loglog(elements, dat[shape]['errors'], label=shape, linestyle='--', marker='x')
ax.set_xlabel('Element Count')
ax.set_ylabel('Relative Error in Solid Body Volume')
ax.set_ylim([0.001, 10.0])
format_spec = lambda l: "{:.3f}\%".format(l) if l < 1.0 else "{:.1f}\%".format(l)
ax.set_yticklabels([format_spec(l) for l in ax.get_yticks().tolist()])
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join('outputs', approach, 'convergence_{:s}.pdf'.format(mesh_type)))
