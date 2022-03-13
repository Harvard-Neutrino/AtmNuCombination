import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time

import util
from params import *
import classdef as cl 



sim = cl.Simulation(pd.read_csv(input_file))
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
fluxes = bf_fluxes
analysis = cl.Analysis(sim, bf_fluxes, fluxes)
util.get_all_weights(analysis, cl.PointType.BestFit)
util.get_all_weights(analysis, cl.PointType.Physical)

numevents = 0
numnuevents = 0
for i in range(2):
	for j in range(2):
		for k in range(2):
			for l in range(len(analysis.bf_weights[i][j][k])):
				if sim.E_tr[l] >= 1 and sim.E_tr[l] <= 50 and np.cos(sim.C_tr[l]) > 0:
					numevents += analysis.bf_weights[i][j][k][l]
					if sim.pdg[l] > 0:
						numnuevents += analysis.bf_weights[i][j][k][l]

print("total number of upgoing events is ", numevents)
print("upgoing neutrinos is ", numnuevents)