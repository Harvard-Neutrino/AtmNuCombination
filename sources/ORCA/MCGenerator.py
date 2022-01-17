import numpy as np
import digitalizer as dgt
import util

class Generator:
	def __init__(self, input_mc_file, input_energy_gaus_track, input_energy_gaus_cas, ORCA_bins):
		self.MC = input_mc_file
		self.G_E_tr = input_energy_gaus_track
		self.G_E_ca = input_energy_gaus_cas
		self.bins = ORCA_bins

	def generate(self, N_IC, N_ORCA):

		# first define function to find reco energy for one event
		def find_reco_energy(gaus, true_energy):
			# set bounds on IC energy
			if true_energy > 53 or true_energy < 1.85:
				return -1

			# find which E_true_ORCA bin this E_true belongs to
			for i in range(len(self.bins) - 1):
				bin_num = -1
				if true_energy <= self.bins[i + 1]:
					bin_num = i
					break

			if bin_num == -1:
				return -1

			# now generate a random reco number from t
			sigma, mu, A = 0, 0, 1
			[sigma, mu, A] = gaus[bin_num]
			random_E_reco = np.random.normal(sigma, mu)

			# return this fake ORCA MC energy
			return random_E_reco

		# now define function to find reco zenith
		def find_reco_zenith(exp, true_energy):
			# bounds on IC energy
			if true_energy > 53 or true_energy < 1.85:
				return -1

			# find sigma
			sigma = exp(true_energy)
			# random generate reco error
			reco_zenith_error = np.random.normal(sigma, 0) * np.pi / 180

			return reco_zenith_error

		# define the 4 exponential functions of zenith error
		e_error, eb_error, mu_error, mub_error = util.get_zenith_error()

		all_e_true = []
		all_e_reco = []
		all_zen_true = []
		all_zen_reco = []
		# now generate a fake ORCA MC energy and ORCA MC weight for all the events
		for i in range(len(self.MC["true_energy"])):
			energy = self.MC["true_energy"][i]
			zenith = self.MC["true_zenith"][i]
			#if i < 10: # just to test the code
			if self.MC["pid"][i] == 0:
				# use cascade gaussian params
				ORCA_E_reco = find_reco_energy(self.G_E_ca, energy)
				ORCA_W_MC = N_ORCA / N_IC * self.MC["weight"][i]
				if int(self.MC["pdg"][i]) > 0:
					ORCA_zen_reco = self.MC["reco_zenith"][i] - find_reco_zenith(e_error, energy)
				else:
					ORCA_zen_reco = self.MC["reco_zenith"][i] - find_reco_zenith(eb_error, energy)
			elif self.MC["pid"][i] == 1:
				# use track gaussian params
				ORCA_E_reco = find_reco_energy(self.G_E_tr, energy)
				ORCA_W_MC = N_ORCA / N_IC * self.MC["weight"][i]
				if int(self.MC["pdg"][i]) > 0:
					ORCA_zen_reco = self.MC["reco_zenith"][i] - find_reco_zenith(mu_error, energy)
				else:
					ORCA_zen_reco = self.MC["reco_zenith"][i] - find_reco_zenith(mub_error, energy)
			else:
				print("invalid pid detected")
				exit(1)

			# print("For event ", i, ":")
			# print("the true energy is ", self.MC["true_energy"][i])
			# print("the IC reco energy is ", self.MC["reco_energy"][i])
			# print("the fake ORCA MC energy is ", ORCA_E_reco, "\n")
			# print("the true zenith is ", self.MC["true_zenith"][i])
			# print("the IC reco zenith is ", self.MC["reco_zenith"][i])
			# print("the fake ORCA MC zenith is ", ORCA_zen_reco, "\n")

			all_e_true.append(energy)
			all_e_reco.append(ORCA_E_reco)
			all_zen_true.append(zenith)
			all_zen_reco.append(ORCA_zen_reco)

		res_e_true = np.array(all_e_true)
		res_e_reco = np.array(all_e_reco)
		res_zen_true = np.array(all_zen_true)
		res_zen_reco = np.array(all_zen_reco)

		return res_e_true, res_e_reco, res_zen_true, res_zen_reco

