import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time
from matplotlib import pyplot as plt

import util
from params import *
import classdef as cl 

def IC_num_events():
    # first check number of total events for 3 years for IC
    sim = cl.Simulation(pd.read_csv(input_file), ORCA = False)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0
    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    num += analysis.bf_weights[i][j][k][l]

    print("IceCube upgrade predicted events per year is ", num / 3)
    return

def ORCA_num_events():
    # first check number of total events for 3 years for IC
    sim = cl.Simulation(pd.read_csv(input_file), ORCA = True)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0
    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    num += analysis.bf_weights[i][j][k][l]

    print("ORCA upgrade predicted events per year is ", num / 3)
    return

IC_num_events()
ORCA_num_events()