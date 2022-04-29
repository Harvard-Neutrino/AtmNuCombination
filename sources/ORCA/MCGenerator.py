import numpy as np
import digitalizer as dgt
from Effective import ICEffectiveAnalysis as ICEff
from Effective import ORCAEffectiveAnalysis as ORCAEff
import util
from scipy.stats import truncnorm
import sys
import math
from params import *

class Generator:
	def __init__(self, input_mc_file, input_energy_gaus_track, input_energy_gaus_cas, ORCA_bins):
		self.MC = input_mc_file
		self.G_E_tr = input_energy_gaus_track
		self.G_E_ca = input_energy_gaus_cas
		self.bins = ORCA_bins

	def generate(self):
		# first define function to find reco energy for one event
		def find_reco_energy(top, true_energy):
			# specify gaussian of topologies
			if top == 0:
				gaus = self.G_E_ca
			elif top == 1:
				gaus = self.G_E_tr
			elif top == 2: # if intermediate
				rand = np.random.randint(pess_ereco)
				if rand == 0:
					gaus = self.G_E_ca
				elif rand >= 1: # this is a pessimistic result
					gaus = self.G_E_tr
			# set bounds on IC energy
			if true_energy >= 54 or true_energy <= 1.85:
				return -1

			# find which E_true_ORCA bin this E_true belongs to
			for i in range(len(self.bins) - 1):
				bin_num = -1
				if true_energy <= self.bins[i + 1]:
					bin_num = i
					break

			if bin_num == -1:
				return -1

			# truncated normal if it belongs to the latter bins (look into truncate function later)
			# now generate a random reco number from t
			sigma, mu, A = 0, 0, 1
			[sigma, mu, A] = gaus[bin_num]
			sigma = np.abs(sigma)
			mu = np.abs(mu)

			# hardcodes the truncated normal distributions, but turned off with lognorm implemented
			# if (top == 1 or top == 2) and bin_num >= 12 and bin_num < 15: # manually set truncate
			# 	while True:
			# 		random_E_reco = np.random.lognormal(sigma, mu)
			# 		if random_E_reco >= 4:
			# 			break
			# elif (top == 1 or top == 2) and bin_num >= 15:
			# 	while True:
			# 		random_E_reco = np.random.lognormal(sigma, mu)
			# 		if random_E_reco >= 5:
			# 			break

			random_E_reco = np.random.lognormal(sigma, mu)

			# return this fake ORCA MC energy
			return random_E_reco

		# now define function to find reco zenith
		def find_reco_zenith(exp, true_energy):
			# bounds on IC energy
			true_energy = min(54, true_energy)
			true_energy = max(1.85, true_energy)

			# find sigma
			sigma = exp(true_energy)
			# random generate reco error
			reco_zenith_error = np.random.normal(0, sigma) * np.pi / 180

			return reco_zenith_error

		# set up effective area + volume analysis
		ICVol = ICEff(1, 50, 51)
		ICVol.computeArea()
		ICVol.computeVolume()
		ORCAVol = ORCAEff()
		ORCAVol.computeArea()
		ORCAVol.computeVolume()
		ICe, ICeb, ICmu, ICmub, ICtau, ICtaub, ICnc, ICncb = ICVol.returnVol()
		ORe, OReb, ORmu, ORmub, ORtau, ORtaub, ORnc, ORncb = ORCAVol.returnVol()


		def find_weight_ratio(true_energy, pdg, interaction_type): # it's actually current type

			volbins = np.logspace(np.log10(1), np.log10(50), num=51)
			
			idx = 0
			for i in range(len(volbins) - 1):
				idx = i
				if true_energy <= volbins[i + 1]:
					break
			
			# print(idx)
			
			if interaction_type == 0:
				if pdg / np.abs(pdg) == 1:
					return ORnc[idx] / ICnc[idx]
				elif pdg / np.abs(pdg) == -1:
					return ORncb[idx] / ICncb[idx]
				else:
					print("wrong pdg detected")
					exit(0)
			elif interaction_type == 1:
				if pdg == 12:
					return ORe[idx] / ICe[idx]
				elif pdg == -12:
					return OReb[idx] / ICeb[idx]
				elif pdg == 14:
					return ORmu[idx] / ICmu[idx]
				elif pdg == -14:
					return ORmub[idx] / ICmub[idx]
				elif pdg == 16:
					return ORtau[idx] / ICtau[idx]
				elif pdg == -16:
					return ORtaub[idx] / ICtaub[idx]
				else:
					print("wrong pdg detected")
					exit(0)
			else:
				print("wrong interaction type detected: ", interaction_type)
				exit(0)

		def assign_topology(nutype, current_type, pdg, true_energy):
			track_prob, cas_prob = util.get_topology_prob(nutype, current_type, pdg, true_energy)

			# now generate discrete choices with numpy
			topologies = [0, 1, 2] # 0 is cascade, 1 is track, 2 is intermediate
			probabilities = [cas_prob, track_prob, 1 - cas_prob - track_prob]
			rand_topology = np.random.choice(topologies, p = probabilities)

			return rand_topology

		# define the 4 exponential functions of zenith error
		e_error, eb_error, mu_error, mub_error = util.get_zenith_error()

		all_e_true = []
		all_e_reco = []
		all_zen_true = []
		all_zen_reco = []
		all_Wmc = []
		all_pid = []
		all_mask = []
		# now generate a fake ORCA MC energy and ORCA MC weight for all the events
		for i in range(len(self.MC["true_energy"])):
			n = len(self.MC["true_energy"]) / 1000
			j = math.floor(i / 1000)
			k = (j + 1) / n
			# print the progress bar
			if i % 1000 == 0:
				print("\r[%-50s] %d%%" % ('='*int(50*k), 100*k), end = '\r', flush = True)

			if self.MC["true_energy"][i] >= 54 or self.MC["true_energy"][i] <= 1.85:
				all_mask.append(False)
				continue

			# first we should actually generate the new topogy, then the energy based on this topology
			# assign random pid's!
			if not rand_morph:
				ORCA_pid = self.MC["pid"][i]
			elif rand_morph:
				ORCA_pid = assign_topology(self.MC["pdg"][i] / np.abs(self.MC["pdg"][i]), self.MC["current_type"][i], \
										self.MC["pdg"][i], self.MC["true_energy"][i]) 

			energy = self.MC["true_energy"][i]
			zenith = self.MC["true_zenith"][i]
			#if i < 10: # just to test the code
			if ORCA_pid == 0:
				# use cascade gaussian params
				if not rand_energy:
					ORCA_E_reco = self.MC["reco_energy"][i]
				elif rand_energy:
					ORCA_E_reco = find_reco_energy(0, energy)
				if not rand_zen:
					ORCA_zen_reco = self.MC["reco_zenith"][i]
				elif rand_zen:
					if int(self.MC["pdg"][i]) > 0:
						ORCA_zen_reco = self.MC["true_zenith"][i] - find_reco_zenith(e_error, energy)
					else:
						ORCA_zen_reco = self.MC["true_zenith"][i] - find_reco_zenith(eb_error, energy)
			elif ORCA_pid == 1:
				# use track gaussian params
				if not rand_energy:
					ORCA_E_reco = self.MC["reco_energy"][i]
				elif rand_energy:
					ORCA_E_reco = find_reco_energy(1, energy)
				if not rand_zen:
					ORCA_zen_reco = self.MC["reco_zenith"][i]
				elif rand_zen:
					if int(self.MC["pdg"][i]) > 0:
						ORCA_zen_reco = self.MC["true_zenith"][i] - find_reco_zenith(mu_error, energy)
					else:
						ORCA_zen_reco = self.MC["true_zenith"][i] - find_reco_zenith(mub_error, energy)
			elif ORCA_pid == 2:
				# use intermediate gaussian params and average on the zenith 
				if not rand_energy:
					ORCA_E_reco = self.MC["reco_energy"][i]
				elif rand_energy:
					ORCA_E_reco = find_reco_energy(2, energy)
				if not rand_zen:
					ORCA_zen_reco = self.MC["reco_zenith"][i]
				elif rand_zen:
					if int(self.MC["pdg"][i]) > 0:
						ORCA_zen_reco = self.MC["true_zenith"][i] - 0.5 * (find_reco_zenith(mu_error, energy) + find_reco_zenith(e_error, energy))
					else:
						ORCA_zen_reco = self.MC["true_zenith"][i] - 0.5 * (find_reco_zenith(mub_error, energy) + find_reco_zenith(eb_error, energy))
			else:
				print("invalid pid detected")
				exit(1)
			
			# also give some new MC weights!
			W_mc = self.MC["weight"][i]
			if not reweight:
				all_Wmc.append(W_mc * 66 / 61)
			elif reweight:
				if energy <= 54 and energy >= 1.85:
					ratio = find_weight_ratio(self.MC["true_energy"][i], self.MC["pdg"][i], self.MC["current_type"][i])
					all_Wmc.append(W_mc * ratio)
				else:
					all_Wmc.append(0)

			all_e_true.append(energy)
			all_e_reco.append(ORCA_E_reco)
			all_zen_true.append(zenith)
			all_zen_reco.append(ORCA_zen_reco)
			all_pid.append(ORCA_pid)
			all_mask.append(True)

		res_e_true = np.array(all_e_true)
		res_e_reco = np.array(all_e_reco)
		res_zen_true = np.array(all_zen_true)
		res_zen_reco = np.array(all_zen_reco)
		res_wmc = np.array(all_Wmc)
		res_pid = np.array(all_pid)

		return res_e_true, res_e_reco, res_zen_true, res_zen_reco, res_wmc, res_pid, all_mask


