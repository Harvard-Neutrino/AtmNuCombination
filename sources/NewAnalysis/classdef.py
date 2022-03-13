import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import sys
from enum import Enum
from scipy.optimize import minimize
 
from params import *

if includeORCA:
	ORCAMC = pd.read_csv("ORCA.csv")

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

class PointType(Enum):
	BestFit = 0
	Physical = 1

class Systematics:
	def __init__(self, norm, nu_nubar_ratio, energy_slope, e_mu_ratio, direction_up, direction_down):
		self.N = norm
		self.delta = nu_nubar_ratio
		self.gamma = energy_slope
		self.eps = e_mu_ratio
		self.Up = direction_up
		self.Down = direction_down

class Simulation:
	def __init__(self, input_file):
		self.file = input_file
		self.W_mc = input_file["weight"]
		self.E_tr = input_file["true_energy"]
		self.C_tr = input_file["true_zenith"]
		self.E_re = input_file["reco_energy"]
		self.C_re = input_file["reco_zenith"]
		self.pdg = input_file["pdg"]
		self.pid = input_file["pid"]
		if includeORCA:
			self.W_mc = ORCAMC["weight"]
			self.E_re = ORCAMC["reco_energy"]
			self.C_re = ORCAMC["reco_zenith"]


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
		self.nsq_atm.Set_CPPhase(0, 2, dcp)

		self.nsq_atm.Set_rel_error(1.0e-4)
		self.nsq_atm.Set_abs_error(1.0e-4)

		self.nsq_atm.Set_initial_state(self.AtmInitFlux,nsq.Basis.flavor)
		self.nsq_atm.EvolveState()

class Analysis:
	def __init__(self, simulation, bf_fluxes, fluxes):
		self.simulation = simulation
		self.bf_fluxes = bf_fluxes
		self.fluxes = fluxes
		# the rest have structures: (neutype 2 x flavor 2 x topology 2)
		self.bf_weights = np.zeros((2, 2, 2, len(self.simulation.W_mc)), dtype = float)
		self.weights = np.zeros((2, 2, 2, len(self.simulation.W_mc)), dtype = float)
		self.bf_histogram = np.zeros((2, 2, 2, NErec, Ncrec), dtype = float)
		self.histogram = np.zeros((2, 2, 2, NErec, Ncrec), dtype = float)
		self.chisq = 0
		self.chisq_min = None
	
	# none-bf part written on the plane and need testing
	def get_weights(self, flavor, neutype, pointtype):
		# print("W_mc: ", self.simulation.W_mc[:10])
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

		# conciseness or performance? if inside or outside the loop
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

			if pointtype.value == 0:
				if self.simulation.E_re[i] > Erec_min and self.simulation.E_tr[i] * units.GeV > E_min and self.simulation.E_tr[i] * units.GeV < E_max:
					self.bf_weights[whatneutype][whatflavor][topology][i] += self.simulation.W_mc[i] * self.bf_fluxes[whatflavor][whatneutype].nsq_atm.EvalFlavor(neuflavor, \
						np.cos(self.simulation.C_tr[i]), self.simulation.E_tr[i] * units.GeV, neutype) * Time * meter_to_cm_sq
			elif pointtype.value == 1:
				if self.simulation.E_re[i] > Erec_min and self.simulation.E_tr[i] * units.GeV > E_min and self.simulation.E_tr[i] * units.GeV < E_max:
					self.weights[whatneutype][whatflavor][topology][i] += self.simulation.W_mc[i] * self.fluxes[whatflavor][whatneutype].nsq_atm.EvalFlavor(neuflavor, \
						np.cos(self.simulation.C_tr[i]), self.simulation.E_tr[i] * units.GeV, neutype) * Time * meter_to_cm_sq
			else:
				print("invalid point type selected")
				exit(1)
		
	# The followings are written on a plane and needs testing


	# This function applies only tilt and zenith
	def pre_apply_systematics(self, sys):

		# in the first step before chi-squared we must first apply slope and cos zenith
		# first the energy slope
		tilt = (self.simulation.E_tr / E0) ** sys.gamma
		# # now the cosine zenith
		cosZen = np.cos(self.simulation.C_tr)
		TanCos = np.tanh(cosZen)
		mask_cUP = cosZen < 0
		mask_cDOWN = cosZen > 0
		zenith = np.ones(len(self.simulation.E_tr))
		zenith[self.mask_cUP] -= sys.Up * self.TanCos[self.mask_cUP]
		zenith[self.mask_cDOWN] -= sys.Down * self.TanCos[self.mask_cDOWN]

		for i in range(2):
			for j in range(2):
				for k in range(2):
					self.weights[i][j][k] *= tilt
					self.weights[i][j][k] *= zenith

	# now we can histogram the weights
	def binning(self, pointtype):
		# performance or conciseness?
		for i in range(2):
			for j in range(2):
				for k in range(2):
					if pointtype.value == 0:
						self.bf_histogram[i][j][k], _, _ = np.histogram2d(self.simulation.E_re, np.cos(self.simulation.C_re), \
													bins = (erec, crec), weights = self.bf_weights[i][j][k])
					elif pointtype.value == 1:
						self.histogram[i][j][k], _, _ = np.histogram2d(self.simulation.E_re, np.cos(self.simulation.C_re), \
													bins = (erec, crec), weights = self.weights[i][j][k])
					else:
						print("invalid point type selected")
						exit(1)

	# now we can calculate chi squared, this part is not completed yet, but it's just the computation of chisq
	def get_chisq(self, sys):

		N = sys.N
		delta = sys.delta 
		eps = sys.eps 

		e_bf = self.bf_histogram[0][0]
		ebar_bf = self.bf_histogram[1][0]
		mu_bf = self.bf_histogram[0][1]
		mubar_bf = self.bf_histogram[1][1]

		e = self.histogram[0][0]
		ebar = self.histogram[1][0]
		mu = self.histogram[0][1]
		mubar = self.histogram[1][1]

		chisq = 0

		for top in range(2):
			for i in range(self.bf_histogram[0][0][0].shape[0]):
				for j in range(self.bf_histogram[0][0][0].shape[1]):
					if e_bf[top][i][j] + ebar_bf[top][i][j] + mu_bf[top][i][j] + mubar_bf[top][i][j] > 0:
						chisq += (N * (delta * eps * e[top][i][j] + (2 - delta) * eps * ebar[top][i][j] + delta * (2 - eps) * mu[top][i][j] + \
																						(2 - delta) * (2 - eps) * mubar[top][i][j]) - \
								(e_bf[top][i][j] + ebar_bf[top][i][j] + mu_bf[top][i][j] + mubar_bf[top][i][j])) ** 2 \
								/ (e_bf[top][i][j] + ebar_bf[top][i][j] + mu_bf[top][i][j] + mubar_bf[top][i][j])

		Sys_BF = Systematics(N_bf, delta_bf, gamma_bf, eps_bf, hv_bf, hv_bf)
		Sys_Sigma = Systematics(sig_N, sig_delta, sig_gamma, sig_eps, sig_hv, sig_hv)

		chisq += (sys.N - Sys_BF.N) ** 2 / Sys_Sigma.N ** 2
		chisq += (sys.delta - Sys_BF.delta) ** 2 / Sys_Sigma.delta ** 2
		chisq += (sys.gamma - Sys_BF.gamma) ** 2 / Sys_Sigma.gamma ** 2
		chisq += (sys.eps - Sys_BF.eps) ** 2 / Sys_Sigma.eps ** 2
		chisq += (sys.Up - Sys_BF.Up) ** 2 / Sys_Sigma.Up ** 2
		chisq += (sys.Down - Sys_BF.Down) ** 2 / Sys_Sigma.Down ** 2

		self.chisq = chisq

		return chisq

	def min_chisq(self):
		# define a function to minimize
		def to_min(syst_array):
			syst = Systematics(syst_array[0], syst_array[1], syst_array[2], syst_array[3], syst_array[4], syst_array[5])
			# self.pre_apply_systematics(syst)
			self.binning(PointType.BestFit)
			self.binning(PointType.Physical)
			self.get_chisq(syst)
			return self.chisq
		
		bnds = ((0., 2.), (0., 2.), (-1., 1), (0., 2.),  (-1., 1), (-1., 1))

		res = minimize(to_min, [N_bf, delta_bf, gamma_bf, eps_bf, hv_bf, hv_bf], bounds = bnds, method = 'L-BFGS-B', options={'disp': True})

		self.chisq_min = res.fun

		return res.fun
