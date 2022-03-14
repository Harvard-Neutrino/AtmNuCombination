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

num_upgoing_events = 0
num_downgoing_events = 0
num_upgoing_nu_events = 0
num_downgoing_nu_events = 0
for i in range(2):
	for j in range(2):
		for k in range(2):
			for l in range(len(analysis.bf_weights[i][j][k])):
				if sim.E_tr[l] >= 1 and sim.E_tr[l] <= 53 and np.cos(sim.C_tr[l]) < 0:
					num_upgoing_events += analysis.bf_weights[i][j][k][l]
					if sim.pdg[l] > 0:
						num_upgoing_nu_events += analysis.bf_weights[i][j][k][l]
				elif sim.E_tr[l] >= 1 and sim.E_tr[l] <= 53 and np.cos(sim.C_tr[l]) > 0:
					num_downgoing_events += analysis.bf_weights[i][j][k][l]
					if sim.pdg[l] > 0:
						num_downgoing_nu_events += analysis.bf_weights[i][j][k][l]

print("total number of upgoing events is ", num_upgoing_events)
print("total number of upgoing neutrinos is ", num_upgoing_nu_events)
print("total number of upgoing events is ", num_downgoing_events)
print("total number of upgoing neutrinos is ", num_downgoing_nu_events)