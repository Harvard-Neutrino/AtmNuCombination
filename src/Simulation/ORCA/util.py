import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, mu, sigma, A):
    return (A / (sigma * np.sqrt(2 * np.pi)) * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2)))

def gaus_fit(data_entries, bins, current_binnum):
	bins_centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
	popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[current_binnum, 5, 1])
	return popt[0], popt[1], popt[2]

def get_zenith_error():

	# first define a general exponential function
	def exp(a, b, c, d, x):
		return a + b * np.exp(c * x) + d * x

	e_neutrino = lambda x : exp(8.29046, 30.2972, -0.3048, -0.12256, x)
	e_antineutrino = lambda x : exp(6.69236, 28.38647, -0.36536, -0.093334, x)
	mu_neutrino = lambda x : exp(8.355, 47.171, -0.45966, -0.10707, x)
	mu_antineutrino = lambda x : exp(6.17314, 42.50309, -0.41, -0.08031, x)

	return e_neutrino, e_antineutrino, mu_neutrino, mu_antineutrino