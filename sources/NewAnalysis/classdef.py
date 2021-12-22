import numpy as np 
import pandas as pd
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux
import sys
from enum import Enum
from scipy.optimize import minimize

from params import *

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

# These still need working on
# define functions to extract lists of enum objects from lists
# def flavor_list(ls):
# 	for i in range(len(ls)):
# 		ls[i] = Flavor(np.absolute(ls[i]))
# 	return ls

# def neutype_list(ls):
# 	for i in range(len(ls)):
# 		ls[i] = NeuType(np.absolute(ls[i])/ls[i])
# 	return ls

# def topology_list(ls):
# 	for i in range(len(ls)):
# 		ls[i] = Topology(ls[i])
# 	return ls 

class Simulation:
	def __init__(self, input_file):
		self.file = input_file
		self.W_mc = input_file["weight"]
		self.E_tr = input_file["true_energy"]
		self.C_tr = input_file["true_zenith"]
		self.E_re = input_file["reco_energy"]
		self.C_re = input_file["reco_zenith"]
		# These need working on
		# self.flavor = flavor_list(input_file["pdg"])
		# self.neutype = neutype_list(input_file["pdg"])
		# self.topology = topology_list(input_file["pid"])
		self.pdg = input_file["pdg"]
		self.pid = input_file["pid"]

class Flux:
	def __init__(self, cth_nodes, energy_nodes):
		self.nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,3,nsq.NeutrinoType.both,interactions)
		self.flux = nuflux.makeFlux("IPhonda2014_spl_solmin")
		self.AtmInitFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))

	def set_initial_flux(self, flavor, neutype):
		def match_flux(flavor, neutype):
			if flavor.value == 12:
				if neutype.value == -1:
					return 1, 0, nuflux.NuEBar
				else: return 0, 0, nuflux.NuE
			elif flavor.value == 14:
				if neutype.value == -1:
					return 1, 1, nuflux.NuMuBar
				else: return 0, 1, nuflux.NuMu
			else:
				print("invalid Atm initial flux selected")
				exit(1)

		what_type, what_flavor, what_flux = match_flux(flavor, neutype)
		for ic,cth in enumerate(self.nsq_atm.GetCosthRange()):
		    for ie,E in enumerate(self.nsq_atm.GetERange()):
		        nu_energy = E/units.GeV
		        nu_cos_zenith = cth
		        self.AtmInitFlux[ic][ie][what_type][what_flavor] = \
		        		self.flux.getFlux(what_flux,nu_energy,nu_cos_zenith)

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
		# self.nsq_atm.Set_ProgressBar(True)
		self.nsq_atm.EvolveState()

class Analysis:
	def __init__(self, simulation, bf_fluxes, fluxes):
		self.simulation = simulation
		self.bf_fluxes = bf_fluxes
		self.fluxes = fluxes
		# self.W_r = np.zeros_like(self.simulation.W_mc)
		self.bf_weights = np.zeros((2, 2, 2, len(self.simulation.W_mc)))
		self.weights = np.zeros((2, 2, 2, len(self.simulation.W_mc)))
	
	def get_bf_weights(self, flavor, neutype):
		def get_flavor_neutype(flavor, neutype):
			if flavor.value == 12:
				if neutype.value == -1:
					return 1, 0
				else: return 0, 0
			elif flavor.value == 14:
				if neutype.value == -1:
					return 1, 1
				else: return 0, 1
			else:
				print("invalid Atm initial flux selected")
				exit(1)
		whatflavor, whatneutype = get_flavor_neutype(flavor, neutype)
	
		for i in range(len(self.simulation.W_mc)):

			if self.simulation.pdg[i] > 0 :
				neutype = 0
			else:
				neutype = 1
						
			if self.simulation.pid[i] == 0 :
				topology = 0
			elif self.simulation.pid[i] == 1 :
				topology = 1

			if np.abs(self.simulation.pdg[i]) == 12:
				neuflavor = 0
			elif np.abs(self.simulation.pdg[i]) == 14:
				neuflavor = 1
			elif np.abs(self.simulation.pdg[i]) == 16:
				neuflavor = 2

			if self.simulation.E_re[i] > Erec_min and self.simulation.E_tr[i] * units.GeV > E_min and self.simulation.E_tr[i] * units.GeV < E_max:
				self.bf_weights[whatneutype][whatflavor][topology][i] += self.simulation.W_mc[i] * self.bf_fluxes[whatflavor][whatneutype].nsq_atm.EvalFlavor(neuflavor, np.cos(self.simulation.C_tr[i]), self.simulation.E_tr[i] * units.GeV, neutype) * Time * meter_to_cm_sq







