import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, mu, sigma):
    return (1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2)))

def gaus_fit(data_entries, bins):
	bins_centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
	# print(bins_centers)
	# print(data_entries)
	popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[0.5, 0.3])
	# print(popt)
	return popt[0], popt[1]

	# # Generate enough x values to make the curves look smooth.
	# xspace = np.linspace(2, 55, 100000)

	# # Plot the histogram and the fitted function.
	# plt.bar(bins_centers, data_entries, width=bins[1] - bins[0], color='navy', label=r'Histogram entries')
	# plt.plot(xspace, gaussian(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')
	# plt.show()