import numpy as np
import scipy as scp
import pandas as pd
from matplotlib import pyplot as plt
from scipy import interpolate

import NewEffective as eff


plt.style.use('./fluxplot.mplstyle')

# first the function to read yerr bar data
def read_yerr(input_file):
	data = pd.read_csv(input_file, header = None, usecols = [0, 1], sep = ',')
	x_val = []
	y_val = []
	yerr_lo = []
	yerr_hi = []
	for i in range(len(data[0])):
		if (i % 3 == 0):
			x_val.append(10**(float(data[0][i])))
			y_val.append(float(data[1][i]))
		elif (i % 3 == 1):
			yerr_hi.append(float(data[1][i] - data[1][i - 1]))
		else:
			yerr_lo.append(float(data[1][i - 2] - data[1][i]))
	yerr = np.array([np.array(yerr_lo),np.array(yerr_hi)])
	return x_val, y_val, yerr

# then the function to read both yerr and xerr bar data
def read_xyerr(input_file):
	data = pd.read_csv(input_file, header = None, usecols = [0, 1], sep = ',')
	x_val = []
	y_val = []
	yerr_lo = []
	yerr_hi = []
	xerr_lo = []
	xerr_hi = []
	for i in range(len(data[0])):
		if (i % 5 == 0):
			x_val.append(10**(float(data[0][i])))
			y_val.append(float(data[1][i]))
		elif (i % 5 == 1):
			yerr_hi.append(float(data[1][i] - data[1][i - 1]))
		elif (i % 5 == 2):
			yerr_lo.append(float(data[1][i - 2] - data[1][i]))
		elif (i % 5 == 3):
			xerr_lo.append(10**(float(data[0][i - 3])) - 10**(float(data[0][i])))
		else:
			xerr_hi.append(10**(float(data[0][i])) - 10**(float(data[0][i - 4])))
	yerr = np.array([np.array(yerr_lo),np.array(yerr_hi)])
	xerr = np.array([np.array(xerr_lo),np.array(xerr_hi)])
	return x_val, y_val, xerr, yerr

# Reads the flux lines as well as interpolates them
def readline(input_file):
	data = pd.read_csv(input_file, header = None, usecols = [0, 1], sep = ',')
	x_val = 10**(data[0])
	y_val = data[1]
	f = interpolate.interp1d(x_val, y_val, kind = "quadratic")
	return f

# Now read in all the flux information
SKnue_x_val, SKnue_y_val, SKnue_yerr = read_yerr("../FluxRangePlot/SKnue.csv")
SKnumu_x_val, SKnumu_y_val, SKnumu_yerr = read_yerr("../FluxRangePlot/SKnumu.csv")
Frej_nue_x_val, Frej_nue_y_val, Frej_nue_xerr, Frej_nue_yerr = read_xyerr("../FluxRangePlot/Frejusnue.csv")
ICDC2013_x_val, ICDC2013_y_val, ICDC2013_xerr, ICDC2013_yerr = read_xyerr("../FluxRangePlot/IC_DC2013nue.csv")
IC2014_x_val, IC2014_y_val, IC2014_xerr, IC2014_yerr = read_xyerr("../FluxRangePlot/IC_DC2013nue.csv")
ANTARES_x_val, ANTARES_y_val, ANTARES_xerr, ANTARES_yerr = read_xyerr("../FluxRangePlot/ANTARES_numu.csv")
Frej_numu_x_val, Frej_numu_y_val, Frej_numu_xerr, Frej_numu_yerr = read_xyerr("../FluxRangePlot/Frejus_numu.csv")
AMANDA_numu_x_val, AMANDA_numu_y_val, AMANDA_numu_xerr, AMANDA_numu_yerr = read_xyerr("../FluxRangePlot/AMANDA-II_numu.csv")
IC_numu_unfolding_x_val, IC_numu_unfolding_y_val, IC_numu_unfolding_xerr, IC_numu_unfolding_yerr = read_xyerr("../FluxRangePlot/IC_numu_unfolding.csv")
HKKM11_nue_f = readline("../FluxRangePlot/HKKM11_nue.csv")
HKKM11_nue_xval = np.logspace(-1, 3.94, 2000)
HKKM11_nue_yval = HKKM11_nue_f(HKKM11_nue_xval)
HKKM11_numu_f = readline("../FluxRangePlot/HKKM11_numu.csv")
HKKM11_numu_xval = np.logspace(-1, 4, 2000)
HKKM11_numu_yval = HKKM11_numu_f(HKKM11_numu_xval)

# Now read in the effective volume information
IC = eff.ICEffectiveAnalysis(0.4, 1.5 * 10 ** 2, 11)
ORCA = eff.ORCAEffectiveAnalysis(num_bins = 10)
IC.compute()
ORCA.compute()
IChist, ICbin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[0][0] + IC.volumes[0][1] + IC.volumes[0][2] \
															+ IC.volumes[1][0]+ IC.volumes[1][1]+ IC.volumes[1][2])
# Now create SK effective volume information
SKbins = np.array([10**-1, 10**(0.7)])
SKweights = np.array([27000])


# start plotting
fig, [ax2, ax] = plt.subplots(figsize = (8, 17), ncols = 1, nrows = 2, gridspec_kw = {'height_ratios':[1, 2]}, sharex = True)
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel(r"Energy [GeV]")
ax.set_ylabel(r"$E^2 \Phi$ [GeV $\cdot$ cm$^{-2}$ sec$^{-1}$ sr$^{-1}]$")
ax2.set_ylabel(r"Effective Volume [m $^3$]")
# ax.set_xlim(-1, 6)
# First plot the flux points
ax.errorbar(SKnue_x_val, SKnue_y_val, yerr = SKnue_yerr, fmt = 'o', capsize = 3, label = r"SK I-IV $\nu_e$")
ax.errorbar(SKnumu_x_val, SKnumu_y_val, yerr = SKnumu_yerr, fmt = 's', capsize = 3, label = r"SK I-IV $\nu_\mu$")
ax.errorbar(Frej_nue_x_val, Frej_nue_y_val, xerr = Frej_nue_xerr, yerr = Frej_nue_yerr, fmt = 'x', capsize = 3, label = r"Frejus $\nu_e$")
ax.errorbar(ICDC2013_x_val, ICDC2013_y_val, xerr = ICDC2013_xerr, yerr = ICDC2013_yerr, fmt = 'v', capsize = 3, label = r"IceCube/Deepcore 2013 $\nu_e$")
ax.errorbar(IC2014_x_val, IC2014_y_val, xerr = IC2014_xerr, yerr = IC2014_yerr, fmt = 'v', capsize = 3, label = r"IceCube 2014 $\nu_e$")
ax.errorbar(ANTARES_x_val, ANTARES_y_val, xerr = ANTARES_xerr, yerr = ANTARES_yerr, fmt = 'P', fillstyle = 'none', capsize = 3, label = r"ANTARES $\nu_mu$")
ax.errorbar(Frej_numu_x_val, Frej_numu_y_val, xerr = Frej_numu_xerr, yerr = Frej_numu_yerr, fmt = 'x', capsize = 3, label = r"Frejus $\nu_\mu$")
ax.errorbar(AMANDA_numu_x_val, AMANDA_numu_y_val, xerr = AMANDA_numu_xerr, yerr = AMANDA_numu_yerr, fmt = 'D', fillstyle = 'none', capsize = 3, label = r"AMANDA-II $\nu_\mu$")
ax.errorbar(IC_numu_unfolding_x_val, IC_numu_unfolding_y_val, xerr = IC_numu_unfolding_xerr, yerr = IC_numu_unfolding_yerr, fmt = '^', capsize = 3, label = r"IceCube $\nu_\mu$ unfolding")
# now plot the fitted flux lines
ax.plot(HKKM11_nue_xval, HKKM11_nue_yval, linestyle = '--', linewidth = 2.4,  label = r"HKKM11 $\nu_e$")
ax.plot(HKKM11_numu_xval, HKKM11_numu_yval, linestyle = '--', linewidth = 2.4, label = r"HKKM11 $\nu_\mu$")

# now label the experiment energy sensitivity regions
ax.axvspan(10**(-1), 10**(0.7), alpha = 0.2)
ax.axvspan(1, 50, alpha = 0.3, color = "lightgreen")
ax.axvspan(0.4, 1.5 * 10 ** 2, alpha = 0.2, color = "orange")

# plt.text(10**(-0.93), 10**(-8), r"Super-K I-IV", fontsize = 18, color = 'steelblue')
# plt.text(10**(-0.93), 10**(-8.3), r" + Super-K Gd", fontsize = 18, color = 'steelblue')

# plt.text(4.5, 10**(-8.1), r"ORCA", fontsize = 18, color = 'seagreen')

# plt.text(1.3, 10**(-8.6), r"IceCube Upgrade", fontsize = 18, color = 'peru')

# Now plot the upper panel
ax2.hist(ICbin[:-1], ICbin, weights = IChist / IC.widths, label=r"IceCube Upgrade Effective Volume", lw = 0, color = "orange", alpha = 0.4)
ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.evol + ORCA.ebarvol + ORCA.muvol + ORCA.mubarvol\
					+ ORCA.tauvol + ORCA.taubarvol + ORCA.ncvol + ORCA.ncbarvol, label=r"ORCA Effective Volume", lw = 0 , color = "lightgreen", alpha = 0.4)
ax2.hist(SKbins[:-1], SKbins, weights = SKweights, label=r"Super-K Effective Volume", lw = 0, alpha = 0.3)

ax2.set_yscale('log')
ax2.legend(fontsize = 14)

# set the x axis range
ax.set_xlim(10**(-1), 1.5 * 10 ** 2)

ax.legend(loc = 3, fontsize = 14)
plt.subplots_adjust(hspace=0.03)
# plt.show()
plt.savefig("./new_paper_plots/Flux_range_plot", bbox_inches = 'tight', pad_inches = 0.4)





