import numpy as np
import digitalizer as dgt

class Generator:
	def __init__(self, input_mc_file, input_energy_gaus_track, input_energy_gaus_cas, ORCA_bins):
		self.MC = input_mc_file
		self.G_E_tr = input_energy_gaus_track
		self.G_E_ca = input_energy_gaus_cas
		self.bins = ORCA_bins

	def generate(self, N_IC, N_ORCA):

		# first define function to find reco energy for one event
		def find_reco(gaus, true_energy):
			# set bounds on IC energy
			if true_energy > 53 or true_energy < 1.85:
				return -1
			# find which E_true_ORCA bin this E_true belongs to
			for i in range(len(self.bins) - 1):
				if true_energy <= self.bins[i + 1]:
					bin_num = i
					break

			# now generate a random reco number from t
			[sigma, mu] = gaus[bin_num]
			random_E_reco = np.random.normal(sigma, mu)

			# return this fake ORCA MC energy
			return random_E_reco

		# now generate a fake ORCA MC energy and ORCA MC weight for all the events
		for i in range(len(self.MC["true_energy"])):
			if i < 10: # just to test the code
				if self.MC["pid"][i] == 0:
					# use cascade gaussian params
					ORCA_E_reco = find_reco(self.G_E_ca, self.MC["true_energy"][i])
					ORCA_W_MC = N_ORCA / N_IC * self.MC["weight"][i]
				elif self.MC["pid"][i] == 1:
					# use track gaussian params
					ORCA_E_reco = find_reco(self.G_E_tr, self.MC["true_energy"][i])
					ORCA_W_MC = N_ORCA / N_IC * self.MC["weight"][i]
				else:
					print("invalid pid detected")
					exit(1)

				print("the true energy is ", self.MC["true_energy"][i])
				print("the fake ORCA MC energy is ", ORCA_E_reco)

