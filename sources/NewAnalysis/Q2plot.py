import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time
from matplotlib import pyplot as plt
import matplotlib

import util
from params import *
import classdef as cl 
# matplotlib.rcParams.update({'font.size': 20})
plt.style.use('./paper.mplstyle')


sim = cl.Simulation(pd.read_csv("neutrino_mc.csv"))
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
fluxes = bf_fluxes

analysis = cl.Analysis(sim, bf_fluxes, fluxes)
util.get_all_weights(analysis, cl.PointType.BestFit)
util.get_all_weights(analysis, cl.PointType.Physical)

cas_weights = np.zeros_like(sim.W_mc)
track_weights = np.zeros_like(sim.W_mc)
for i in range(len(cas_weights)):
    cas_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][1][0][i] +  \
                        analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][1][0][i]
    track_weights[i] = analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][1][1][i] +  \
                        analysis.bf_weights[1][0][1][i] + analysis.bf_weights[1][1][1][i]

Q2bins = np.logspace(-4.7, 1, 21)
bin_widths = np.zeros((20,))
for i in range(20):
    bin_widths[i] = Q2bins[i+1] - Q2bins[i]

cascade, _ = np.histogram(sim.file["Q2"], bins = Q2bins, weights = cas_weights)
track, _ = np.histogram(sim.file["Q2"], bins = Q2bins, weights = track_weights)

fig, ax = plt.subplots(figsize = (8, 8))

ax.hist(Q2bins[:-1], Q2bins, weights = cascade / bin_widths,\
                    label="Cascade", histtype="step")
ax.hist(Q2bins[:-1], Q2bins, weights = track / bin_widths,\
                    label="Track", histtype="step")

ax.legend()
ax.set_xlabel(r"$Q^2$")
ax.set_ylabel("Event Count")
ax.set_xscale("log")
ax.set_yscale("log")
plt.show()
fig.savefig("ICQ2plot")