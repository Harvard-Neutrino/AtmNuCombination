import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import pandas as pd
from numpy import sin 
from numpy import cos  


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


# define the rodrigues rotation formula
def rand_vector():
	ph = np.random.uniform(low = 0, high = np.pi * 2)
	th = np.random.uniform(low = 0, high = np.pi)
	# print(ph, th)
	u = np.array([cos(ph)*sin(th),sin(ph)*sin(th),cos(th)])
	# print(u)
	# u = np.array([1, 1, 0])
	u = u / np.linalg.norm(u)
	return u

# def the conversion between coordinate systems
def cart_from_sp(th, ph):
	u = np.array([cos(ph)*sin(th),sin(ph)*sin(th),cos(th)])
	return u

# directly get the rotation matrix
def R_from_axis(u, th):
	cth = cos(th)
	sth = sin(th)
	R11 = cos(th) + u[0]**2*(1-cos(th))
	R12 = u[0]*u[1]*(1-cos(th)) - u[2]*sin(th)
	R13 = u[0] * u[2] * (1-cos(th)) + u[1] * sin(th)
	R21 = u[1] * u[0] * (1-cos(th)) + u[2] * sth
	R22 = cth + u[1]**2 * (1-cth)
	R23 = u[1] * u[2] * (1-cth) - u[0] * sth
	R31 = u[2] * u[0] * (1-cth) - u[1] * sth
	R32 = u[2] * u[1] * (1-cth) + u[0] * sth
	R33 = cth + u[2] ** 2 * (1-cth)
	return np.array([[R11, R12, R13],[R21, R22, R23],[R31, R32, R33]])

# checks the theta is retained by trace
def check_trace(R):
	Tr = R[0][0] + R[1][1] + R[2][2]
	return np.arccos((Tr - 1) / 2)


# implementation using rodrigues rotation formula
def rod_rot(v, u, ang):
	u = u / np.linalg.norm(u)
	vrot = np.array([0.,0.,0.])
	vrot[0] = (cos(ang)+u[0]**2*(1-cos(ang)))*v[0] + (u[0]*u[1]*(1-cos(ang))-u[2]*sin(ang))*v[1] + (u[0]*u[2]*(1-cos(ang)+u[1]*sin(ang)))*v[2]
	vrot[1] = (cos(ang)+u[1]**2*(1-cos(ang)))*v[1] + (u[0]*u[1]*(1-cos(ang))+u[2]*sin(ang))*v[0] + (u[1]*u[2]*(1-cos(ang)-u[0]*sin(ang)))*v[2]
	vrot[2] = (cos(ang)+u[2]**2*(1-cos(ang)))*v[2] + (u[0]*u[2]*(1-cos(ang))-u[1]*sin(ang))*v[0] + (u[1]*u[2]*(1-cos(ang)+u[0]*sin(ang)))*v[1]
	# print(vrot)
	vrot = vrot / np.linalg.norm(vrot)
	print(v[0] * vrot[0] + v[1] * vrot[1] + v[2] * vrot[2])
	# print(vrot[0] ** 2 + vrot[1] ** 2 + vrot[2] ** 2) # checked that it's normalized
	return vrot

# implementation using the two rotations
def rot_vector(v, th, phi):
	v = v / np.linalg.norm(v)
	# first create an orthogonal vector
	u = np.array([v[1], -v[0], 0])
	u = u / np.linalg.norm(u)

	# rotate u around this vector by theta
	Rth = R_from_axis(u, th)

	vth = np.zeros((3,))
	vth[0] = Rth[0][0] * v[0] + Rth[0][1] * v[1] + Rth[0][2] * v[2]
	vth[1] = Rth[1][0] * v[0] + Rth[1][1] * v[1] + Rth[1][2] * v[2]
	vth[2] = Rth[2][0] * v[0] + Rth[2][1] * v[1] + Rth[2][2] * v[2]
	vth = vth / np.linalg.norm(vth)

	# now apply the phi rotation in a cone
	Rphi = R_from_axis(v, phi)
	vrot = np.zeros((3,))
	vrot[0] = Rphi[0][0] * vth[0] + Rphi[0][1] * vth[1] + Rphi[0][2] * vth[2]
	vrot[1] = Rphi[1][0] * vth[0] + Rphi[1][1] * vth[1] + Rphi[1][2] * vth[2]
	vrot[2] = Rphi[2][0] * vth[0] + Rphi[2][1] * vth[1] + Rphi[2][2] * vth[2]
	vrot = vrot / np.linalg.norm(vrot)
	return vrot

# v = np.array([1, 0, 0])
# print(rot_vector(v, np.pi / 2))


# check the dot product
def check_dot(v, th, phi):
	vrot = rot_vector(v, th, phi)
	return np.arccos(np.dot(v, vrot)/(np.linalg.norm(v) * np.linalg.norm(vrot)))

# print(check_dot(rand_vector(), 2, np.random.uniform(low = 0, high = 2 * np.pi)))

# this version is the rodrigues rotation formula implementation
# def rand_reco_zen(zen, azim, err):
# 	v = np.array([np.cos(azim) * np.sin(zen), np.sin(azim) * np.sin(zen), np.cos(zen)])
# 	# print(v)
# 	u = rand_vector()
# 	vrot = rod_rot(v, u, err)
# 	zenrot = (np.arctan(np.sqrt(vrot[0] ** 2 + vrot[1] ** 2) / vrot[2]))
# 	if zenrot < 0:
# 		zenrot += np.pi
# 	return zenrot

# this is the two successive rotations version
def rand_reco_zen(zen, azim, err):
	v = np.array([np.cos(azim) * np.sin(zen), np.sin(azim) * np.sin(zen), np.cos(zen)])
	# print(v)
	vrot = rot_vector(v, err, np.random.uniform(low = 0, high = 2 * np.pi))
	zenrot = (np.arctan(np.sqrt(vrot[0] ** 2 + vrot[1] ** 2) / vrot[2]))
	if zenrot < 0:
		zenrot += np.pi
	return zenrot


def test_rot_zen_distribution(num):
	ls = np.zeros(num)
	for i in range(num):
		ls[i] = rand_reco_zen(1.5, 2, np.pi/6)

	plt.hist(ls, bins = 30)
	plt.show()
	plt.close()

# test_rot_zen_distribution(100)

def test_rot_graphic(zen, azim, err, num):
	ls = np.ndarray((num, 3))
	zenls = np.zeros(num)
	v = cart_from_sp(zen, azim)
	for i in range(num):
		vrot = rot_vector(v, err, np.random.uniform(low = 0, high = 2 * np.pi))
		ls[i][0] = vrot[0]
		ls[i][1] = vrot[1]
		ls[i][2] = vrot[2]
		zenls[i] = (np.arctan(np.sqrt(vrot[0] ** 2 + vrot[1] ** 2) / vrot[2]))
		if zenls[i] < 0:
			zenls[i] += np.pi
	x_pos = y_pos = z_pos = np.zeros(num)
	x_dir = np.zeros(num)
	y_dir = np.zeros(num)
	z_dir = np.zeros(num)
	for i in range(num):
		x_dir[i] = ls[i][0]
		y_dir[i] = ls[i][1]
		z_dir[i] = ls[i][2]
	fig = plt.figure()
	ax = plt.axes(projection = '3d')
	ax.set_xlim(-1, 1)
	ax.set_ylim(-1, 1)
	ax.set_zlim(-1, 1)
	# start = [0, 0, 0]
	# for i in range(num):
	# 	ax.quiver(start[0], start[1], start[2], ls[i][0], ls[i][1], ls[i][2])
	ax.quiver(x_pos, y_pos, z_pos, x_dir, y_dir, z_dir)
	plt.show()
	plt.close()
	plt.hist(zenls, bins = 30)
	plt.show()
	plt.close()

# test_rot_graphic(1.5, 2, np.pi / 6, 100)



def ORCA_paper_sensitivity():
	data = pd.read_csv("./ORCA_Results/ORCA_sensitivity.csv", header = None, usecols = [0, 1])
	# print(data[0])
	# print(data[1])
	x = np.array(data[0])
	y = np.array(data[1])
	# print(x)
	# fit splines to x=f(u) and y=g(u), treating both as periodic. also note that s=0
	# is needed in order to force the spline fit to pass through all the input points.
	tck, u = interpolate.splprep([x, y], s=0, k=2, per=True)

	# evaluate the spline fits for 1000 evenly spaced distance values
	xi, yi = interpolate.splev(np.linspace(0, 1, 1000), tck)

	# plot the result
	fig, ax = plt.subplots(1, 1)
	# ax.plot(x, y, 'or')
	ax.plot(xi, yi, '-b')

	# plt.show()
	plt.savefig("./ORCA_Results/ORCA_paper_sensitivity")


# ORCA_paper_sensitivity()






















































