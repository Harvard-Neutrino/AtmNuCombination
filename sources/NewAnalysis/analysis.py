import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux

import util
from params import *
import classdef as cl 

# define a simulation from a csv file
sim = cl.Simulation(pd.read_csv(input_file))

# initialize and propagate all 4 = 2x2 fluxes
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)

# in this trial we do fluxes = bf_fluxes
fluxes = bf_fluxes

# start the analysis
analysis = cl.Analysis(sim, bf_fluxes, fluxes)

# get all the rated weight of 2x2 neutrino types
util.get_all_weights(analysis, cl.PointType.BestFit)
util.get_all_weights(analysis, cl.PointType.Physical)

# apply zenith and tilt
analysis.pre_apply_systematics(Syst)

# bin them
analysis.histogram(cl.PointType.BestFit)
analysis.histogram(cl.PointType.Physical)

# apply the other systematics and calculate chi squared
analysis.get_chisq(Syst)