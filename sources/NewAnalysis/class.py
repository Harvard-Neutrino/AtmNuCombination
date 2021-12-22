import numpy as np 
import pandas as pd
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux
import sys
from enum import Enum
from scipy.optimize import minimize

import params as pm 

class Topology(Enum):
    cascade = 0
    track = 1

class NeuType(Enum):
    AntiNeutrino = -1
    Neutrino = 1

class Flavor(Enum):
    e = 12
    mu = 14
    tau = 16

# define functions to extract lists of enum objects from lists
def flavor_list(ls):
	for i in range(len(ls)):
		ls[i] = Flavor(np.absolute(ls[i]))
	return ls

def neutype_list(ls):
	for i in range(len(ls)):
		ls[i] = NeuType(np.absolute(ls[i])/ls[i])
	return ls

def topology_list(ls):
	for i in range(len(ls)):
		ls[i] = Topology(ls[i])
	return ls 

class Simulation:
	def __init__(self, input_file):
		self.file = input_file
		self.W_mc = input_file["weight"]
		self.E_tr = input_file["true_energy"]
		self.C_tr = input_file["true_zenith"]
		self.E_re = input_file["reco_energy"]
		self.C_re = input_file["reco_zenith"]
		self.flavor = flasvor_list(input_file["pdg"])
		self.neutype = neutype_list(input_file["pdg"])
		self.topology = topology_list(input_file["pid"])

class Flux:
	def __init__(self, cth_nodes, energy_nodes):
		self.nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,3,nsq.NeutrinoType.both,interactions)
		self.nuflux = nuflux.makeflux("IPhonda2014_spl_solmin")
		self.AtmInitFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))

	def set_initial_flux(self, flavor, neutype):
		def match_flux(flavor, neutype):
			match flavor.value:
				case 12:
					if neutype.value == -1
						return 1, 0, nuflux.NuEBar
					else: return 0, 0, nuflux.NuE
				case 14:
					if neutype.value == -1:
						return 1, 1, nuflux.NuMuBar
					else: return 0, 1, nuflux.NuMu
				case _:
					print("invalid Atm initial flux selected")
					exit(1)

		what_type, what_flavor, what_flux = match_flux(flavor, neutype)
		for ic,cth in enumerate(nsq_atm.GetCosthRange()):
		    for ie,E in enumerate(nsq_atm.GetERange()):
		        nu_energy = E/units.GeV
		        nu_cos_zenith = cth
		        self.AtmInitFlux[ic][ie][what_type][what_flavor] = \
		        		flux.getFlux(what_flux,nu_energy,nu_cos_zenith)

	def propagate_flux(self, t23, dm31, dcp):
		self.nsq_atm.Set_MixingAngle(0, 1, 33.44 * Radian)
		self.nsq_atm.Set_MixingAngle(0, 2, 8.57 * Radian)
		self.nsq_atm.Set_MixingAngle(1, 2, t23)
		self.nsq_atm.Set_SquareMassDifference(1, 7.42e-5)
		self.nsq_atm.Set_SquareMassDifference(2, dm31)
		self.nsq_atm.Set_CPPhase(0,2, dcp)

		self.nsq_atm.Set_rel_error(1.0e-4)
		self.nsq_atm.Set_abs_error(1.0e-4)

		self.nsq_atm.Set_initial_state(self.AtmInitFlux,nsq.Basis.flavor)
		self.nsq_atm.EvolveState()

class Analysis:
	def __init__(self, simulation, bf_fluxes, fluxes):
		self.simulation = simulation
		self.bf_fluxes = bf_fluxes
		self.fluxes = fluxes
		self.W_r = np.zeros_like(self.simulation.W_mc)
		self.weights = np.zeros(2, len(self.W_r), 2, 2)





