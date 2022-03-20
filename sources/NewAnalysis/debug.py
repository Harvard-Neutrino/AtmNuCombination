import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time

import util
from params import *
import classdef as cl 



sim = cl.Simulation(pd.read_csv(input_file))
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
fluxes = bf_fluxes
analysis = cl.Analysis(sim, bf_fluxes, fluxes)
util.get_all_weights(analysis, cl.PointType.BestFit)
util.get_all_weights(analysis, cl.PointType.Physical)


def print_ORCA_num_event_details():
	num_upgoing_events = 0
	num_downgoing_events = 0
	num_upgoing_nu_events = 0
	num_downgoing_nu_events = 0
	for i in range(2):
		for j in range(2):
			for k in range(2):
				for l in range(len(analysis.bf_weights[i][j][k])):
					if sim.E_tr[l] >= 1 and sim.E_tr[l] <= 53 and np.cos(sim.C_tr[l]) < 0:
						num_upgoing_events += analysis.bf_weights[i][j][k][l]
						if sim.pdg[l] > 0:
							num_upgoing_nu_events += analysis.bf_weights[i][j][k][l]
					elif sim.E_tr[l] >= 1 and sim.E_tr[l] <= 53 and np.cos(sim.C_tr[l]) > 0:
						num_downgoing_events += analysis.bf_weights[i][j][k][l]
						if sim.pdg[l] > 0:
							num_downgoing_nu_events += analysis.bf_weights[i][j][k][l]

	print("total number of upgoing events is ", num_upgoing_events)
	print("total number of upgoing neutrinos is ", num_upgoing_nu_events)
	print("total number of upgoing events is ", num_downgoing_events)
	print("total number of upgoing neutrinos is ", num_downgoing_nu_events)
	return

def ORCA_event_topology_details():
	# define the bins
	energy_bins = np.logspace(0, np.log10(50), 31)
	# first do the numu neutrino tracks
	numu_CC_tracks_weights = np.zeros(len(analysis.bf_weights[0][0][1]))
	numu_CC_all_weights = np.zeros(len(analysis.bf_weights[0][0][1]))
	for i in range(len(numu_tracks_weights)):
		numu_CC_all_weights[i] = analysis.bf_weights[0][0][0][i] + \
					analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][0][2][i]
		if sim.currents[i] == 1: # charged current
			# choose all the numu_CC events
			numu_CC_tracks_weights[i] = analysis.bf_weights[0][0][1][i]

	# now do the fraction histograms
	numu_CC_track_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numu_CC_tracks_weights)
	numu_CC_all_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numu_CC_all_weights)

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("ORCA MC nu_mu CC topology fractions")
	ax.hist(energy_bins[:-1], energy_bins, weights = numu_CC_track_hist / numu_CC_all_hist,\
					 label=r"$\nu_{\mu}CC$", color="blue", histtype="step")

	ax.set_xscale("log")
	ax.set_xlabel("neutrino energy [GeV]")
	ax.set_xlim(1, 50)
	ax.set_ylabel("fraction")
	ax.grid(True)
	ax.legend(loc = 2)
	plt.show()
	# fig.savefig("./RecoPlots/IC_Effective_volume")
	plt.close()

ORCA_event_topology_details()
