import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time

import util
from params import *
import classdef as cl 

startprogram = time.time()

Syst = cl.Systematics(N_bf, delta_bf, gamma_bf, eps_bf, hv_bf, hv_bf)

# define a simulation from a csv file
sim = cl.Simulation(pd.read_csv(input_file))

# initialize and propagate all 4 = 2x2 fluxes
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)

# in this trial we do fluxes = bf_fluxes
fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, np.arcsin(np.sqrt(0.500)), dm31, dcp)

finishprop = time.time()

# start the analysis
analysis = cl.Analysis(sim, bf_fluxes, fluxes)

# get all the rated weight of 2x2 neutrino types
util.get_all_weights(analysis, cl.PointType.Physical)
util.get_all_weights(analysis, cl.PointType.BestFit)

numevents = 0
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(len(analysis.simulation.W_mc)):
                if analysis.simulation.E_tr[l] <= 53 and analysis.simulation.E_tr[l] >= 1.85:
                    numevents += analysis.bf_weights[i][j][k][l]

print(numevents / 3)
exit(0)


finishweights = time.time()

# analysis.prepare_syst()

print(analysis.min_chisq())

finishmin = time.time()

print("propagation took ", finishprop - startprogram)
print("calculation of rated weight took ", finishweights - finishprop)
print("minimization took ", finishmin - finishweights)
print("entire process took ", finishmin - startprogram)
'''

# The followings are for calculating a single chi-squared value from an input syst array

# apply zenith and tilt
analysis.pre_apply_systematics(Syst)

finishpreapply = time.time()

# bin them
analysis.binning(cl.PointType.BestFit)
analysis.binning(cl.PointType.Physical)

finishbinning = time.time()


# apply the other systematics and calculate chi squared
print(analysis.get_chisq(Syst))

finishprogram = time.time()

print("propagation took ", finishprop - startprogram)
print("calculation of rated weight took ", finishweights - finishprop)
print("preapply systematics took ", finishpreapply - finishweights)
print("binning took ", finishbinning - finishpreapply)
print("rest of syst and chisq took ", finishprogram - finishbinning)
print("entire program took ", finishprogram - startprogram)

'''