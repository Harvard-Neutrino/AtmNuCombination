import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time
from matplotlib import pyplot as plt
import matplotlib

import util
from params import *
import classdef as cl 
matplotlib.rcParams.update({'font.size': 20})



sim = cl.Simulation(pd.read_csv("../ORCA/15x_with_interm.csv"))
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
fluxes = bf_fluxes

analysis = cl.Analysis(sim, bf_fluxes, fluxes)
util.get_all_weights(analysis, cl.PointType.BestFit)
util.get_all_weights(analysis, cl.PointType.Physical)


# plot the energy distribution of tracks and cascades
def true_event_distribution():
	cas_weights = np.zeros_like(sim.W_mc)
	track_weights = np.zeros_like(sim.W_mc)
	for i in range(len(cas_weights)):
		cas_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][1][0][i] +  \
							analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][1][0][i]
		track_weights[i] = analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][1][1][i] +  \
							analysis.bf_weights[1][0][1][i] + analysis.bf_weights[1][1][1][i]
	energy_bins = np.logspace(0, np.log10(50), 21)
	bin_widths = np.zeros((20,))
	for i in range(20):
		bin_widths[i] = energy_bins[i+1] - energy_bins[i]
	cascade, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = cas_weights)
	track, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = track_weights)


	zen_bins = np.linspace(-1, 0, 11)
	zen_bin_widths = np.zeros((10,))
	for i in range(10):
		zen_bin_widths[i] = zen_bins[i+1] - zen_bins[i]
	# print(bin_widths)
	zencascade, _ = np.histogram(np.cos(sim.C_tr), bins = zen_bins, weights = cas_weights)
	zentrack, _ = np.histogram(np.cos(sim.C_tr), bins = zen_bins, weights = track_weights)
	# print(cascade)



	fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(20,10))
	ax1, ax2 = axes[0], axes[1]
	fig.suptitle("ORCA MC Event Distributions")

	ax2.hist(zen_bins[:-1], zen_bins, weights = zencascade / zen_bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax2.hist(zen_bins[:-1], zen_bins, weights = zentrack / zen_bin_widths,\
					 label="ORCA track", histtype="step")
	# ax.set_xscale("log")
	ax2.set_yscale("log")
	ax2.set_xlabel("Neutrino True Cosine Zenith")
	ax2.set_xlim(-1, 0)
	ax2.set_ylabel("Total Event Count [3yrs]")
	# ax.grid(True)
	ax2.legend()
	ax2.title.set_text("True Cosine Zenith Distribution")
	ax1.title.set_text("True Energy Distribution")

	ax1.hist(energy_bins[:-1], energy_bins, weights = cascade / bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax1.hist(energy_bins[:-1], energy_bins, weights = track / bin_widths,\
					 label="ORCA track", histtype="step")
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.set_xlabel("Neutrino True Energy [GeV]")
	ax1.set_xlim(1.85, 50)
	ax1.set_ylabel("Total Event Count [3yrs]")
	ax1.set_xticks([1.85, 10, 50])
	# ax.grid(True)
	ax1.legend()
	# plt.show()
	fig.savefig("./../ORCA/new_paper_plots/Event_Distribution_True.png")
	plt.close()

# true_event_distribution()


# plot the reco distribution of tracks and cascades
def reco_event_distribution():
	cas_weights = np.zeros_like(sim.W_mc)
	track_weights = np.zeros_like(sim.W_mc)
	for i in range(len(cas_weights)):
		cas_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][1][0][i] +  \
							analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][1][0][i]
		track_weights[i] = analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][1][1][i] +  \
							analysis.bf_weights[1][0][1][i] + analysis.bf_weights[1][1][1][i]
	energy_bins = np.logspace(0, np.log10(50), 21)
	bin_widths = np.zeros((20,))
	for i in range(20):
		bin_widths[i] = energy_bins[i+1] - energy_bins[i]
	cascade, _ = np.histogram(sim.E_re, bins = energy_bins, weights = cas_weights)
	track, _ = np.histogram(sim.E_re, bins = energy_bins, weights = track_weights)


	zen_bins = np.linspace(-1, 0, 11)
	zen_bin_widths = np.zeros((10,))
	for i in range(10):
		zen_bin_widths[i] = zen_bins[i+1] - zen_bins[i]
	# print(bin_widths)
	zencascade, _ = np.histogram(np.cos(sim.C_re), bins = zen_bins, weights = cas_weights)
	zentrack, _ = np.histogram(np.cos(sim.C_re), bins = zen_bins, weights = track_weights)
	# print(cascade)



	fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(20,10))
	ax1, ax2 = axes[0], axes[1]
	# fig.suptitle("ORCA MC Event Distributions")

	ax2.hist(zen_bins[:-1], zen_bins, weights = zencascade / zen_bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax2.hist(zen_bins[:-1], zen_bins, weights = zentrack / zen_bin_widths,\
					 label="ORCA track", histtype="step")
	# ax.set_xscale("log")
	ax2.set_yscale("log")
	ax2.set_xlabel("Neutrino Reconstructed Cosine Zenith")
	ax2.set_xlim(-1, 0)
	ax2.set_ylabel("Total Event Count [3yrs]")
	# ax.grid(True)
	ax2.legend(loc = 2)
	ax2.title.set_text("Reconstructed Cosine Zenith Distribution")
	ax1.title.set_text("Reconstructed Energy Distribution")

	ax1.hist(energy_bins[:-1], energy_bins, weights = cascade / bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax1.hist(energy_bins[:-1], energy_bins, weights = track / bin_widths,\
					 label="ORCA track", histtype="step")
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.set_ylim(2 * 10 ** 2, 2 * 10 ** 4)
	ax1.set_xlabel("Neutrino Reconstructed Energy [GeV]")
	ax1.set_xlim(1.85, 50)
	ax1.set_ylabel("Total Event Count [3yrs]")
	ax1.set_xticks([1.85, 10, 50])
	# ax.grid(True)
	ax1.legend()
	# plt.show()
	fig.savefig("./../ORCA/new_paper_plots/Event_Distribution_Reco.png")
	plt.close()

reco_event_distribution()

def true_event_distribution_zen():
	cas_weights = np.zeros_like(sim.W_mc)
	track_weights = np.zeros_like(sim.W_mc)
	for i in range(len(cas_weights)):
		cas_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][1][0][i] +  \
							analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][1][0][i]
		track_weights[i] = analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][1][1][i] +  \
							analysis.bf_weights[1][0][1][i] + analysis.bf_weights[1][1][1][i]
	energy_bins = np.logspace(0, np.log10(50), 21)
	bin_widths = np.zeros((20,))
	for i in range(20):
		bin_widths[i] = energy_bins[i+1] - energy_bins[i]
	cascade, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = cas_weights)
	track, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = track_weights)


	zen_bins = np.linspace(-0.3, np.pi + 0.3, 11)
	zen_bin_widths = np.zeros((10,))
	for i in range(10):
		zen_bin_widths[i] = zen_bins[i+1] - zen_bins[i]
	# print(bin_widths)
	zencascade, _ = np.histogram(sim.C_tr, bins = zen_bins, weights = cas_weights)
	zentrack, _ = np.histogram(sim.C_tr, bins = zen_bins, weights = track_weights)
	# print(cascade)



	fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(20,10))
	ax1, ax2 = axes[0], axes[1]
	fig.suptitle("ORCA MC Event Distributions")

	ax2.hist(zen_bins[:-1], zen_bins, weights = zencascade / zen_bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax2.hist(zen_bins[:-1], zen_bins, weights = zentrack / zen_bin_widths,\
					 label="ORCA track", histtype="step")
	# ax.set_xscale("log")
	ax2.set_yscale("log")
	ax2.set_xlabel("Neutrino True Zenith")
	ax2.set_xlim(-0.3, np.pi + 0.3)
	ax2.set_ylabel("Total Event Count [3yrs]")
	# ax.grid(True)
	ax2.legend()
	ax2.title.set_text("True Zenith Distribution")
	ax1.title.set_text("True Energy Distribution")

	ax1.hist(energy_bins[:-1], energy_bins, weights = cascade / bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax1.hist(energy_bins[:-1], energy_bins, weights = track / bin_widths,\
					 label="ORCA track", histtype="step")
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.set_xlabel("Neutrino True Energy [GeV]")
	ax1.set_xlim(1.85, 50)
	ax1.set_ylabel("Total Event Count [3yrs]")
	ax1.set_xticks([1.85, 10, 50])
	# ax.grid(True)
	ax1.legend()
	# plt.show()
	fig.savefig("./../ORCA/new_paper_plots/Event_Distribution_True_w_zenith.png")
	plt.close()

# true_event_distribution_zen()


# plot the reco distribution of tracks and cascades
def reco_event_distribution_zen():
	cas_weights = np.zeros_like(sim.W_mc)
	track_weights = np.zeros_like(sim.W_mc)
	for i in range(len(cas_weights)):
		cas_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][1][0][i] +  \
							analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][1][0][i]
		track_weights[i] = analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][1][1][i] +  \
							analysis.bf_weights[1][0][1][i] + analysis.bf_weights[1][1][1][i]
	energy_bins = np.logspace(0, np.log10(50), 21)
	bin_widths = np.zeros((20,))
	for i in range(20):
		bin_widths[i] = energy_bins[i+1] - energy_bins[i]
	cascade, _ = np.histogram(sim.E_re, bins = energy_bins, weights = cas_weights)
	track, _ = np.histogram(sim.E_re, bins = energy_bins, weights = track_weights)


	zen_bins = np.linspace(-0.3, np.pi + 0.3, 11)
	zen_bin_widths = np.zeros((10,))
	for i in range(10):
		zen_bin_widths[i] = zen_bins[i+1] - zen_bins[i]
	# print(bin_widths)
	zencascade, _ = np.histogram(sim.C_re, bins = zen_bins, weights = cas_weights)
	zentrack, _ = np.histogram(sim.C_re, bins = zen_bins, weights = track_weights)
	# print(cascade)



	fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(20,10))
	ax1, ax2 = axes[0], axes[1]
	fig.suptitle("ORCA MC Event Distributions")

	ax2.hist(zen_bins[:-1], zen_bins, weights = zencascade / zen_bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax2.hist(zen_bins[:-1], zen_bins, weights = zentrack / zen_bin_widths,\
					 label="ORCA track", histtype="step")
	# ax.set_xscale("log")
	ax2.set_yscale("log")
	ax2.set_xlabel("Neutrino Reco Zenith")
	ax2.set_xlim(-0.3, np.pi + 0.3)
	ax2.set_ylabel("Total Event Count [3yrs]")
	# ax.grid(True)
	ax2.legend()
	ax2.title.set_text("Reco Zenith Distribution")
	ax1.title.set_text("Reco Energy Distribution")

	ax1.hist(energy_bins[:-1], energy_bins, weights = cascade / bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax1.hist(energy_bins[:-1], energy_bins, weights = track / bin_widths,\
					 label="ORCA track", histtype="step")
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.set_xlabel("Neutrino Reco Energy [GeV]")
	ax1.set_xlim(1.85, 50)
	ax1.set_ylabel("Total Event Count [3yrs]")
	ax1.set_xticks([1.85, 10, 50])
	# ax.grid(True)
	ax1.legend()
	# plt.show()
	fig.savefig("./../ORCA/new_paper_plots/Event_Distribution_Reco_w_zenith.png")
	plt.close()

# reco_event_distribution_zen()