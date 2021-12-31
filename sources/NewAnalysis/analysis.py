import numpy as np 
import pandas as pd
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux
import sys
from enum import Enum
from scipy.optimize import minimize

from params import *
import classdef as cl 

sim = cl.Simulation(pd.read_csv(input_file))

bf_e_flux = cl.Flux(cth_nodes, energy_nodes)
bf_e_flux.set_initial_flux(cl.Flavor.e, cl.NeuType.Neutrino)
bf_e_flux.propagate_flux(theta23, dm31, dcp)

bf_mu_flux = cl.Flux(cth_nodes, energy_nodes)
bf_mu_flux.set_initial_flux(cl.Flavor.mu, cl.NeuType.Neutrino)
bf_mu_flux.propagate_flux(theta23, dm31, dcp)

bf_ebar_flux = cl.Flux(cth_nodes, energy_nodes)
bf_ebar_flux.set_initial_flux(cl.Flavor.e, cl.NeuType.AntiNeutrino)
bf_ebar_flux.propagate_flux(theta23, dm31, dcp)

bf_mubar_flux = cl.Flux(cth_nodes, energy_nodes)
bf_mubar_flux.set_initial_flux(cl.Flavor.mu, cl.NeuType.AntiNeutrino)
bf_mubar_flux.propagate_flux(theta23, dm31, dcp)

bf_fluxes = np.array([[bf_e_flux, bf_ebar_flux],[bf_mu_flux, bf_mubar_flux]])

# in this trial we do fluxes = bf_fluxes
fluxes = np.array([[bf_e_flux, bf_ebar_flux],[bf_mu_flux, bf_mubar_flux]])

analysis = cl.Analysis(sim, bf_fluxes, fluxes)

# Rest are not tested

# should also include util functions to bundle these functions together

analysis.get_weights()

analysis.pre_apply_systematics()

analysis.histogram()

analysis.get_chisq()