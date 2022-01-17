import sys
import os
import numpy as np

import util
import classdef as cl 
from params import *

# this file calculates a chi-squared value for a pair of physical parameters
# and then saves the physical parameters as well as the chi-sq into a file

idx = int(sys.argv[1])

# get the physical parameters
t23val, m31val = util.id_to_param(idx)

print("t23 and m31 are: {}, {}".format(t23val, m31val))

file_name = "{}.txt".format(idx)
complete_name = os.path.join(dir_name, file_name)

# define a simulation from a csv file
sim = cl.Simulation(pd.read_csv(input_file))

# initialize and propagate all 4 = 2x2 fluxes
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)

# in this trial we do fluxes = bf_fluxes
fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, t23val, m31val, dcp)

# start the analysis
analysis = cl.Analysis(sim, bf_fluxes, fluxes)

# get all the rated weight of 2x2 neutrino types
util.get_all_weights(analysis, cl.PointType.Physical)
util.get_all_weights(analysis, cl.PointType.BestFit)

chisqval = analysis.min_chisq()

saveres = np.array([idx, chisqval])

np.savetxt(complete_name, saveres)