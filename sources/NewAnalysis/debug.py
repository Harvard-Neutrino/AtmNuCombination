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
for i in range(2):
	for j in range(2):
		for k in range(2):
			for l in range(len(analysis.bf_weights[i][j][k])):
			 numevents += analysis.bf_weights[i][j][k][l]

print(numevents)