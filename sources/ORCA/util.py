import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import pandas as pd

# gaussian function
def gaussian(x, mu, sigma, A):
    return (A / (sigma * np.sqrt(2 * np.pi)) * np.exp(-1.0 * (np.log(x) - mu)**2 / (2 * sigma**2)))

# fit with the gaussian function
def gaus_fit(data_entries, bins, current_binnum):
	bins_centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
	popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[current_binnum, 5, 1])
	return popt[0], popt[1], popt[2]

############################ This does not work for the moment#######################################################
# two gaussian function
def twogaussian(x, mu1, mu2, sigma1, sigma2, alpha):
	first = (alpha / (sigma1 * np.sqrt(2 * np.pi)) * np.exp(-1.0 * (np.log(x) - mu1)**2 / (2 * sigma1**2)))
	second = ((1 - alpha)) / (sigma2 * np.sqrt(2 * np.pi)) * np.exp(-1.0 * (np.log(x) - mu2)**2 / (2 * sigma2**2))
	return first + second

# fit with the two gaussian function
def two_gaus_fit(data_entries, bins, current_binnum):
	bins_centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
	popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[current_binnum, 5, 1])
	return popt[0], popt[1], popt[2] , popt[3] , popt[4]
######################################################################################################################


# tanh fit
def tanh(x, A, b, k):
	return b + A * np.tanh(k * (x))

def tanh_fit(data, bins_centers):
	popt, pcov = curve_fit(tanh, xdata=bins_centers, ydata=data, p0=[1.5, 1, 2.8])
	return popt[0], popt[1], popt[2]

# input the zenith error histogram from ORCA paper plot
def get_zenith_error():

	# first define a general exponential function
	def exp(a, b, c, d, x):
		return a + b * np.exp(c * x) + d * x

	e_neutrino = lambda x : exp(8.29046, 30.2972, -0.3048, -0.12256, x)
	e_antineutrino = lambda x : exp(6.69236, 28.38647, -0.36536, -0.093334, x)
	mu_neutrino = lambda x : exp(8.355, 47.171, -0.45966, -0.10707, x)
	mu_antineutrino = lambda x : exp(6.17314, 42.50309, -0.41, -0.08031, x)

	return e_neutrino, e_antineutrino, mu_neutrino, mu_antineutrino

def getORCAbins(input_file, tau = False, nc = False):
	df = pd.read_csv(input_file, header = None, usecols = [1])
	# print(df)
	# this is the bin heights
	w_bin = np.array(df[:]).T[0] * 10 ** 6
	if tau:
		for i in range(len(w_bin)):
			if i <= 15:
				w_bin[i] = 0
	if nc:
		for i in range(len(w_bin)):
			if i <= 3:
				w_bin[i] = 0
	# print(res)
	return w_bin

# def interp_tau_xs(nutype):
# 	tau_xs = pd.read_csv("./ORCA_Results/nutau_xs.csv", header = None, usecols = [0, 1])
# 	# nubarxsection = pd.read_csv("nubar.txt", sep = ' ', usecols = [0, 1, 2])
# 	if nutype == 1:
# 		length = len(tau_xs[0])
# 		extracted = np.zeros(length)
# 		energies = np.zeros(length)
# 		for i in range(length):
# 			extracted[i] = tau_xs[1][i]
# 			energies[i] = tau_xs[0][i]

# 	resf = interp1d(energies, extracted)

# 	return resf

# # interpolates the tau xs ratio file
# def tau_xs_ratio():
# 	tau_ratio = pd.read_csv("./ORCA_Results/tau_xs_ratio.csv", header = None, usecols = [0, 1])

# 	length = len(tau_ratio[0])
# 	extracted = np.zeros(length)
# 	energies = np.zeros(length)
# 	for i in range(length):
# 		extracted[i] = tau_ratio[1][i]
# 		energies[i] = tau_ratio[0][i]

# 	resf = interp1d(energies, extracted)

# 	return resf

def interpolate_xsection(nutype, tau = False):
	# reads the xsection txt file
	nuxsection = pd.read_csv("nu.txt", sep = ' ', usecols = [0, 1, 2])
	nubarxsection = pd.read_csv("nubar.txt", sep = ' ', usecols = [0, 1, 2])
	# tausec = interp_tau_xs(1)
	# tau_ratio = tau_xs_ratio()
	# taubarsec

	length = len(nuxsection["Energy"])
	extracted = np.zeros(length)
	energies = np.zeros(length)

	if nutype == 1:
		for i in range(length):
			extracted[i] = nuxsection["sigmaCC"][i] + nuxsection["sigmaNC"][i]
			energies[i] = nuxsection["Energy"][i]
		# if not tau:
		# 	for i in range(length):
		# 		extracted[i] = nuxsection["sigmaCC"][i] + nuxsection["sigmaNC"][i]
		# 		energies[i] = nuxsection["Energy"][i]
		# 		# print(extracted[i])
		# 		# print(energies[i])
		# elif tau:
		# 	for i in range(length):
		# 		energy = energies[i]
		# 		if energy >= 3.5 and energy <= 70:
		# 			extracted[i] = tausec(energy) + nuxsection["sigmaNC"][i]
		# 			# extracted[i] = nuxsection["sigmaCC"][i] * tau_ratio(energy) + nubarxsection["sigmaNC"][i]
		# 			energies[i] = nuxsection["Energy"][i]
		# 		else:
		# 			extracted[i] = nuxsection["sigmaCC"][i] + nubarxsection["sigmaNC"][i]
		# 			energies[i] = nuxsection["Energy"][i]
	elif nutype == -1:
		for i in range(length):
			extracted[i] = nubarxsection["sigmaCC"][i] + nubarxsection["sigmaNC"][i]
			energies[i] = nubarxsection["Energy"][i]

	resf = interp1d(energies, extracted)

	return resf

def all_xsec():
	xsec = pd.read_csv("all_xsecs.csv")
	energy = xsec["energy"]
	nu = xsec["nu"]
	nubar = xsec["nubar"]
	tau = xsec["tau"]
	taubar = xsec["taubar"]
	nc = xsec["nc"]
	ncbar = xsec["ncbar"]
	f_nu = interp1d(energy, nu)
	f_nubar = interp1d(energy, nubar)
	f_tau = interp1d(energy, tau)
	f_taubar = interp1d(energy, taubar)
	f_nc = interp1d(energy, nc)
	f_ncbar = interp1d(energy, ncbar)
	return f_nu, f_nubar, f_tau, f_taubar, f_nc, f_ncbar


# return the index of effective volume/area matrix according to current type and pdg
# first index is neutrino/antineutrino = 0/1
# second index is flavor+current (e/mu/tau/nc) = 0/1/2/3
def get_index(current, pdg):
	if current == 1: #charged
		if pdg / np.abs(pdg) == 1: # neutrino
			if pdg == 12:
				return 0, 0
			elif pdg == 14:
				return 0, 1
			elif pdg == 16:
				return 0, 2
		elif pdg / np.abs(pdg) == -1: # antineutrino
			if np.abs(pdg) == 12:
				return 1, 0
			elif np.abs(pdg) == 14:
				return 1, 1
			elif np.abs(pdg) == 16:
				return 1, 2
	if current == 0: # neutral
		if pdg / np.abs(pdg) == 1: # neutrino
			return 0, 3
		elif pdg / np.abs(pdg) == -1: # antineutrino
			return 1, 3


def get_ORCA_topology_prob(nutype, current_type, pdg, true_energy):

	# first define how to get the probabilities
	def get_probs(input_file):
		df = pd.read_csv(input_file, header = None, usecols = [1])
		# these are the histogram weights
		res = np.array(df[:]).T[0]
		if len(res) < 30:
			# print("yes")
			padding = np.zeros(30 - len(res))
			res = np.concatenate((padding, res), axis = 0)
		return res

	# return the index in the topology histogram given the energy
	def energy_to_index(energy):
		energies = np.logspace(np.log10(2), np.log10(50), 31)
		idx = 0
		for i in range(len(energies) - 1):
			# print(i)
			# print(energy)
			# print(energies[i+1])
			if energy <= energies[i + 1]:
				return idx
			idx += 1
		return 29

	if current_type == 0: # NC
		if nutype == 1: # neutrinos
			track = get_probs("./ORCA_Results/track_nu_NC.csv")
			cascade = get_probs("./ORCA_Results/cascade_nu_NC.csv")
		elif nutype == -1: # antineutrinos
			track = get_probs("./ORCA_Results/track_nubar_NC.csv")
			cascade = get_probs("./ORCA_Results/cascade_nubar_NC.csv")

	if current_type == 1: # CC
		if np.abs(pdg) == 12: # nu_e
			if nutype == 1: # neutrinos
				track = get_probs("./ORCA_Results/track_nue_CC.csv")
				cascade = get_probs("./ORCA_Results/cascade_nue_CC.csv")
			elif nutype == -1:
				track = get_probs("./ORCA_Results/track_nuebar_CC.csv")
				cascade = get_probs("./ORCA_Results/cascade_nuebar_CC.csv")
		elif np.abs(pdg) == 14: #nu_mu
			if nutype == 1: # neutrinos
				track = get_probs("./ORCA_Results/track_numu_CC.csv")
				cascade = get_probs("./ORCA_Results/cascade_numu_CC.csv")
			elif nutype == -1:
				track = get_probs("./ORCA_Results/track_numubar_CC.csv")
				cascade = get_probs("./ORCA_Results/cascade_numubar_CC.csv")
		elif np.abs(pdg) == 16: # nu_tau
			if nutype == 1: # neutrinos
				track = get_probs("./ORCA_Results/track_nutau_CC.csv")
				cascade = get_probs("./ORCA_Results/cascade_nutau_CC.csv")
			elif nutype == -1:
				track = get_probs("./ORCA_Results/track_nutaubar_CC.csv")
				cascade = get_probs("./ORCA_Results/cascade_nutaubar_CC.csv")
	
	# now use energy to give two numbers from the two arrays
	idx = energy_to_index(true_energy)
	return track[idx], 1 - cascade[idx]

def get_IC_topology_prob(nutype, current_type, pdg, true_energy):
	def energy_to_index(energy):
		energies = np.logspace(np.log10(2), np.log10(50), 31)
		idx = 0
		for i in range(len(energies) - 1):
			# print(i)
			# print(energy)
			# print(energies[i+1])
			if energy <= energies[i + 1]:
				return idx
			idx += 1
		return 29
	
	def get_prob(nutype, filename, padding = False):
		df = pd.read_csv(filename)
		if nutype == 1: # neutrinos
			track = df["nu_track"]
			cascade = df["nu_cas"]
		elif nutype == -1:
			track = df["nubar_track"]
			cascade = df["nubar_cas"]
		else:
			print("wrong nutype detected")
			exit(1)
		if padding:
			pad = np.zeros(30 - len(track))
			track = np.concatenate((padding, track), axis = 0)
			cascade = np.concatenate((padding, cascade), axis = 0)
		return track, cascade
	
	tau_padding = False
	if current_type == 0: # NC
		filename = "./ORCA_Results/nu_NC_Topology_Fraction"
	elif current_type == 1:
		if np.abs(pdg) == 12:
			filename = "./ORCA_Results/nue_CC_Topology_Fraction"
		elif np.abs(pdg) == 14:
			filename = "./ORCA_Results/numu_CC_Topology_Fraction"
		elif np.abs(pdg) == 16:
			filename = "./ORCA_Results/nutau_CC_Topology_Fraction"
			tau_padding = True
	else:
		print("wrong pdg detected")
	
	track, cascade = get_prob(nutype, filename)

	idx = energy_to_index(true_energy)
	return track[idx], cascade[idx]

# print(get_topology_prob(1, 1, 14, 3))
# print(get_IC_topology_prob(1, 0, 14, 3))


