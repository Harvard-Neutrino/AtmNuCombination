import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time
from matplotlib import pyplot as plt

import util
from params import *
import classdef as cl 
 


sim = cl.Simulation(pd.read_csv("../ORCA/15x_with_interm.csv"))
ic_sim = cl.Simulation(pd.read_csv("neutrino_mc.csv"))
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
fluxes = bf_fluxes

analysis = cl.Analysis(sim, bf_fluxes, fluxes)
util.get_all_weights(analysis, cl.PointType.BestFit)
util.get_all_weights(analysis, cl.PointType.Physical)


ic_analysis = cl.Analysis(ic_sim, bf_fluxes, fluxes)
util.get_all_weights(ic_analysis, cl.PointType.BestFit)
util.get_all_weights(ic_analysis, cl.PointType.Physical)

# checge the flavor basis so we now discern by the 3 flavors in detection
def change_flavor_basis(weights):
	res = np.zeros((2, 3, 3, len(sim.W_mc)))
	for i in range(len(sim.W_mc)):
		a = sim.pdg[i] / np.abs(sim.pdg[i])
		if a == 1:
			neutype = 0
		elif a == -1:
			neutype = 1
		topology = int(sim.pid[i])
		if np.abs(sim.pdg[i]) == 12: # if reco as a nu e
			res[neutype][0][topology][i] = weights[neutype][0][topology][i] + weights[neutype][1][topology][i]
		elif np.abs(sim.pdg[i]) == 14: # if reco as a nu mu
			res[neutype][1][topology][i] = weights[neutype][0][topology][i] + weights[neutype][1][topology][i]
		elif np.abs(sim.pdg[i]) == 16: # if reco as a nu tau
			res[neutype][2][topology][i] = weights[neutype][0][topology][i] + weights[neutype][1][topology][i]
	return res
			

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

def ORCA_topology_details(flavor, current):
	legendloc = 1
	# define the names I save it to
	if flavor == 0:
		flavname = "e"
		labelname1 = r"$\nu_e$ CC Track"
		labelname2 = r"$\nu_e$ CC Cascade"
		labelname3 = r"$\bar{\nu}_e$ CC Track"
		labelname4 = r"$\bar{\nu}_e$ CC Cascade"
	elif flavor == 1:
		flavname = "mu"
		labelname1 = r"$\nu_{\mu}$ CC Track"
		labelname2 = r"$\nu_{\mu}$ CC Cascade"
		labelname3 = r"$\bar{\nu}_{\mu}$ CC Track"
		labelname4 = r"$\bar{\nu}_{\mu}$ CC Cascade"
		legendloc = 4
	elif flavor == 2:
		flavname = "tau"
		labelname1 = r"$\nu_{\tau}$ CC Track"
		labelname2 = r"$\nu_{\tau}$ CC Cascade"
		labelname3 = r"$\bar{\nu}_{\tau}$ CC Track"
		labelname4 = r"$\bar{\nu}_{\tau}$ CC Cascade"
	
	if current == 1:
		curname = "CC"
	elif current == 0:
		curname = "NC"
		flavname = ""
		labelname1 = r"$\nu$ NC Track"
		labelname2 = r"$\nu$ NC Cascade"
		labelname3 = r"$\bar{\nu}$ NC Track"
		labelname4 = r"$\bar{\nu}$ NC Cascade"
	# transform to reco flavor basis
	weights = change_flavor_basis(analysis.bf_weights)

	# define the bins
	energy_bins = np.logspace(np.log10(2), np.log10(50), 31)
	if flavor == 2: # if we look at taus
		energy_bins = np.logspace(np.log10(4), np.log10(50), 31)
	one = np.ones(30)
	mclen = len(analysis.bf_weights[0][1][1])
	# first do the numu neutrino tracks
	numu_CC_tracks_weights = np.zeros(mclen)
	numu_CC_cas_weights = np.zeros(mclen)
	numu_CC_inter_weights = np.zeros(mclen)
	numu_CC_all_weights = np.zeros(mclen)
	numubar_CC_tracks_weights = np.zeros(mclen)
	numubar_CC_cas_weights = np.zeros(mclen)
	numubar_CC_inter_weights = np.zeros(mclen)
	numubar_CC_all_weights = np.zeros(mclen)
	numuavg_CC_tracks_weights = np.zeros(mclen)
	numuavg_CC_cas_weights = np.zeros(mclen)
	numuavg_CC_inter_weights = np.zeros(mclen)
	numuavg_CC_all_weights = np.zeros(mclen)
	for i in range(len(numu_CC_tracks_weights)):
		if current == 1:
			if sim.currents[i] == 1: # charged current
				# choose all the numu_CC events
				numu_CC_all_weights[i] = weights[0][flavor][0][i] + \
						weights[0][flavor][1][i] + weights[0][flavor][2][i]
				numu_CC_tracks_weights[i] = weights[0][flavor][1][i]
				numu_CC_cas_weights[i] = weights[0][flavor][0][i]
				numu_CC_inter_weights[i] = weights[0][flavor][2][i]
				# choose all the numubar CC events
				numubar_CC_all_weights[i] = weights[1][flavor][0][i] + \
						weights[1][flavor][1][i] + weights[1][flavor][2][i]
				numubar_CC_tracks_weights[i] = weights[1][flavor][1][i]
				numubar_CC_cas_weights[i] = weights[1][flavor][0][i]
				numubar_CC_inter_weights[i] = weights[1][flavor][2][i]
				# now the average ones
				numuavg_CC_all_weights[i] = numu_CC_all_weights[i] + numubar_CC_all_weights[i]
				numuavg_CC_tracks_weights[i] = numu_CC_tracks_weights[i] + numubar_CC_tracks_weights[i]
				numuavg_CC_cas_weights[i] = numu_CC_cas_weights[i] + numubar_CC_cas_weights[i]
				numuavg_CC_inter_weights[i] = numu_CC_inter_weights[i] + numubar_CC_inter_weights[i]
		elif current == 0:
			if sim.currents[i] == 0: # neutral current
				for flavor in range(3):
					numu_CC_all_weights[i] += weights[0][flavor][0][i] + \
							weights[0][flavor][1][i] + weights[0][flavor][2][i]
					numu_CC_tracks_weights[i] += weights[0][flavor][1][i]
					numu_CC_cas_weights[i] += weights[0][flavor][0][i]
					numu_CC_inter_weights[i] += weights[0][flavor][2][i]
					# choose all the numubar CC events
					numubar_CC_all_weights[i] += weights[1][flavor][0][i] + \
							weights[1][flavor][1][i] + weights[1][flavor][2][i]
					numubar_CC_tracks_weights[i] += weights[1][flavor][1][i]
					numubar_CC_cas_weights[i] += weights[1][flavor][0][i]
					numubar_CC_inter_weights[i] += weights[1][flavor][2][i]
					# now the average ones
					numuavg_CC_all_weights[i] += numu_CC_all_weights[i] + numubar_CC_all_weights[i]
					numuavg_CC_tracks_weights[i] += numu_CC_tracks_weights[i] + numubar_CC_tracks_weights[i]
					numuavg_CC_cas_weights[i] += numu_CC_cas_weights[i] + numubar_CC_cas_weights[i]
					numuavg_CC_inter_weights[i] += numu_CC_inter_weights[i] + numubar_CC_inter_weights[i]
	
	# now do the fraction histograms
	numu_CC_track_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numu_CC_tracks_weights)
	numu_CC_cas_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numu_CC_cas_weights)
	numu_CC_inter_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numu_CC_inter_weights)
	numu_CC_all_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numu_CC_all_weights)

	numubar_CC_track_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numubar_CC_tracks_weights)
	numubar_CC_cas_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numubar_CC_cas_weights)
	numubar_CC_inter_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numubar_CC_inter_weights)
	numubar_CC_all_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numubar_CC_all_weights)
	
	numuavg_CC_track_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numuavg_CC_tracks_weights)
	numuavg_CC_cas_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numuavg_CC_cas_weights)
	numuavg_CC_inter_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numuavg_CC_inter_weights)
	numuavg_CC_all_hist, _ = np.histogram(sim.E_tr, bins = energy_bins, \
							weights = numuavg_CC_all_weights)
	
	# # this is to save the bin heights (for IC, but general purpose)
	# df_array = np.ndarray((4, 30), dtype = float) # in the order of nucas, nutrack, nubarcas, nubartrack
	# df_array[0] = numu_CC_cas_hist / numu_CC_all_hist
	# df_array[1] = numu_CC_track_hist / numu_CC_all_hist
	# df_array[2] = numubar_CC_cas_hist / numubar_CC_all_hist
	# df_array[3] = numubar_CC_track_hist / numubar_CC_all_hist
	# df = pd.DataFrame(df_array.T, columns = ["nu_cas", "nu_track", "nubar_cas", "nubar_track"])
	# df.to_csv("nu{}_{}_Topology_Fraction".format(flavname, curname))

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle(r"ORCA MC $\nu_\{}$ {} topology fractions".format(flavname, curname))
	# first the neutrino ones
	ax.hist(energy_bins[:-1], energy_bins, weights = numu_CC_track_hist / numu_CC_all_hist,\
					 label=labelname1, color="blue", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = one - numu_CC_cas_hist / numu_CC_all_hist,\
					 label=labelname2, color="red", histtype="step")
	
	# now the antineutrino ones
	ax.hist(energy_bins[:-1], energy_bins, weights = numubar_CC_track_hist / numubar_CC_all_hist,\
					 label=labelname3, linestyle = "--", color="blue", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = one - numubar_CC_cas_hist / numubar_CC_all_hist,\
					 label=labelname4, linestyle = "--", color="red", histtype="step")
	
	#now the average ones
	ax.hist(energy_bins[:-1], energy_bins, weights = one, color="lightpink")
	ax.hist(energy_bins[:-1], energy_bins, weights = one - numuavg_CC_cas_hist / numuavg_CC_all_hist, color="lightgray")
	ax.hist(energy_bins[:-1], energy_bins, weights = numuavg_CC_track_hist / numuavg_CC_all_hist, color="cornflowerblue")
	

	ax.set_xscale("log")
	ax.set_xlabel("neutrino energy [GeV]")
	ax.set_xlim(2, 50)
	ax.set_ylim(0, 1)
	ax.set_ylabel("fraction")
	ax.grid(True)
	ax.legend(loc = legendloc)
	# plt.show()
	fig.savefig("./../ORCA/new_paper_plots/nu{}_{}_Topology_Fraction".format(flavname, curname))
	plt.close()

ORCA_topology_details(0, 1)
ORCA_topology_details(1, 1)
ORCA_topology_details(2, 1)
ORCA_topology_details(3, 0)

# plot the energy distribution of tracks and cascades
def morph_distribution():
	cas_weights = np.zeros_like(sim.W_mc)
	track_weights = np.zeros_like(sim.W_mc)
	ic_cas_weights = np.zeros_like(ic_sim.W_mc)
	ic_track_weights = np.zeros_like(ic_sim.W_mc)
	for i in range(len(cas_weights)):
		cas_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][1][0][i] +  \
							analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][1][0][i]
		track_weights[i] = analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][1][1][i] +  \
							analysis.bf_weights[1][0][1][i] + analysis.bf_weights[1][1][1][i]
	for i in range(len(ic_cas_weights)):
		ic_cas_weights[i] = ic_analysis.bf_weights[0][0][0][i] + ic_analysis.bf_weights[0][1][0][i] +  \
							ic_analysis.bf_weights[1][0][0][i] + ic_analysis.bf_weights[1][1][0][i]
		ic_track_weights[i] = ic_analysis.bf_weights[0][0][1][i] + ic_analysis.bf_weights[0][1][1][i] +  \
							ic_analysis.bf_weights[1][0][1][i] + ic_analysis.bf_weights[1][1][1][i]	
	energy_bins = np.logspace(0, np.log10(50), 21)
	bin_widths = np.zeros((20,))
	for i in range(20):
		bin_widths[i] = energy_bins[i+1] - energy_bins[i]
	# print(bin_widths)
	cascade, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = cas_weights)
	track, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = track_weights)
	ic_cascade, _ = np.histogram(ic_sim.E_tr, bins = energy_bins, weights = ic_cas_weights)
	ic_track, _ = np.histogram(ic_sim.E_tr, bins = energy_bins, weights = ic_track_weights)
	# print(cascade)
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("ORCA MC energy distribution")
	ax.hist(energy_bins[:-1], energy_bins, weights = cascade / bin_widths,\
					 label="ORCA cascade", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = track / bin_widths,\
					 label="ORCA track", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = ic_cascade / bin_widths,\
					 label="IC cascade", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = ic_track / bin_widths,\
					 label="IC track", histtype="step")
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel("neutrino energy [GeV]")
	ax.set_xlim(1, 50)
	ax.set_ylabel("event rate [3yrs]")
	# ax.grid(True)
	ax.legend()
	# plt.show()
	fig.savefig("./../ORCA/RecoPlots/0621poly_rw_distribution_morph")
	plt.close()

# morph_distribution()

# plot the energy distribution of tracks and cascades
def flavor_distribution():
	e_weights = np.zeros_like(sim.W_mc)
	mu_weights = np.zeros_like(sim.W_mc)
	tau_weights = np.zeros_like(sim.W_mc)

	ic_e_weights = np.zeros_like(ic_sim.W_mc)
	ic_mu_weights = np.zeros_like(ic_sim.W_mc)
	ic_tau_weights = np.zeros_like(ic_sim.W_mc)

	for i in range(len(e_weights)):
		if np.abs(sim.pdg[i]) == 12:
			e_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][0][1][i] +  \
								analysis.bf_weights[0][1][0][i] + analysis.bf_weights[0][1][1][i] + \
								analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][0][1][i] +  \
								analysis.bf_weights[1][1][0][i] + analysis.bf_weights[1][1][1][i]
		elif np.abs(sim.pdg[i]) == 14: 
			mu_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][0][1][i] +  \
								analysis.bf_weights[0][1][0][i] + analysis.bf_weights[0][1][1][i] + \
								analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][0][1][i] +  \
								analysis.bf_weights[1][1][0][i] + analysis.bf_weights[1][1][1][i]
		elif np.abs(sim.pdg[i]) == 16: 
			tau_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][0][1][i] +  \
								analysis.bf_weights[0][1][0][i] + analysis.bf_weights[0][1][1][i] + \
								analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][0][1][i] +  \
								analysis.bf_weights[1][1][0][i] + analysis.bf_weights[1][1][1][i]
	for i in range(len(ic_e_weights)):
		if np.abs(ic_sim.pdg[i]) == 12:
			ic_e_weights[i] = ic_analysis.bf_weights[0][0][0][i] + ic_analysis.bf_weights[0][0][1][i] +  \
								ic_analysis.bf_weights[0][1][0][i] + ic_analysis.bf_weights[0][1][1][i] + \
								ic_analysis.bf_weights[1][0][0][i] + ic_analysis.bf_weights[1][0][1][i] +  \
								ic_analysis.bf_weights[1][1][0][i] + ic_analysis.bf_weights[1][1][1][i]
		elif np.abs(ic_sim.pdg[i]) == 14: 
			ic_mu_weights[i] = ic_analysis.bf_weights[0][0][0][i] + ic_analysis.bf_weights[0][0][1][i] +  \
								ic_analysis.bf_weights[0][1][0][i] + ic_analysis.bf_weights[0][1][1][i] + \
								ic_analysis.bf_weights[1][0][0][i] + ic_analysis.bf_weights[1][0][1][i] +  \
								ic_analysis.bf_weights[1][1][0][i] + ic_analysis.bf_weights[1][1][1][i] 
		elif np.abs(ic_sim.pdg[i]) == 16: 
			ic_tau_weights[i] = ic_analysis.bf_weights[0][0][0][i] + ic_analysis.bf_weights[0][0][1][i] +  \
								ic_analysis.bf_weights[0][1][0][i] + ic_analysis.bf_weights[0][1][1][i] + \
								ic_analysis.bf_weights[1][0][0][i] + ic_analysis.bf_weights[1][0][1][i] +  \
								ic_analysis.bf_weights[1][1][0][i] + ic_analysis.bf_weights[1][1][1][i]
	energy_bins = np.logspace(0, np.log10(50), 21)
	bin_widths = np.zeros((20,))
	for i in range(20):
		bin_widths[i] = energy_bins[i+1] - energy_bins[i]
	# print(bin_widths)
	e, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = e_weights)
	mu, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = mu_weights)
	tau, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = tau_weights)

	ic_e, _ = np.histogram(ic_sim.E_tr, bins = energy_bins, weights = ic_e_weights)
	ic_mu, _ = np.histogram(ic_sim.E_tr, bins = energy_bins, weights = ic_mu_weights)
	ic_tau, _ = np.histogram(ic_sim.E_tr, bins = energy_bins, weights = ic_tau_weights)
	# print(cascade)
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("ORCA MC energy distribution")
	ax.hist(energy_bins[:-1], energy_bins, weights = e / bin_widths,\
					 label="ORCA nue", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = mu / bin_widths,\
					 label="ORCA numu", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = tau / bin_widths,\
					 label="ORCA nutau", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = ic_e / bin_widths,\
					 label="IC nue", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = ic_mu / bin_widths,\
					 label="IC numu", histtype="step")
	ax.hist(energy_bins[:-1], energy_bins, weights = ic_tau / bin_widths,\
					 label="IC nutau", histtype="step")
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel("neutrino energy [GeV]")
	ax.set_xlim(1, 50)
	ax.set_ylabel("event rate [3yrs]")
	# ax.grid(True)
	ax.legend()
	# plt.show()
	fig.savefig("./../ORCA/RecoPlots/0621poly_rw_distribution_flavor")
	plt.close()

# flavor_distribution()

def y_distribution():
	# cas_weights = np.zeros_like(sim.W_mc)
	# track_weights = np.zeros_like(sim.W_mc)
	ic_cas_weights = np.zeros_like(ic_sim.W_mc)
	ic_track_weights = np.zeros_like(ic_sim.W_mc)
	# for i in range(len(cas_weights)):
	# 	cas_weights[i] = analysis.bf_weights[0][0][0][i] + analysis.bf_weights[0][1][0][i] +  \
	# 						analysis.bf_weights[1][0][0][i] + analysis.bf_weights[1][1][0][i]
	# 	track_weights[i] = analysis.bf_weights[0][0][1][i] + analysis.bf_weights[0][1][1][i] +  \
	# 						analysis.bf_weights[1][0][1][i] + analysis.bf_weights[1][1][1][i]
	for i in range(len(ic_cas_weights)):
		ic_cas_weights[i] = ic_analysis.bf_weights[0][0][0][i] + ic_analysis.bf_weights[0][1][0][i] +  \
							ic_analysis.bf_weights[1][0][0][i] + ic_analysis.bf_weights[1][1][0][i]
		ic_track_weights[i] = ic_analysis.bf_weights[0][0][1][i] + ic_analysis.bf_weights[0][1][1][i] +  \
							ic_analysis.bf_weights[1][0][1][i] + ic_analysis.bf_weights[1][1][1][i]	
	binnum = 10
	y_bins = np.logspace(-4, 0, binnum + 1)
	bin_widths = np.zeros((binnum,))
	for i in range(binnum):
		bin_widths[i] = y_bins[i+1] - y_bins[i]
	# print(bin_widths)
	# cascade, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = cas_weights)
	# track, _ = np.histogram(sim.E_tr, bins = energy_bins, weights = track_weights)
	ic_cascade, _ = np.histogram(ic_sim.file["y"], bins = y_bins, weights = ic_cas_weights)
	ic_track, _ = np.histogram(ic_sim.file["y"], bins = y_bins, weights = ic_track_weights)
	# print(cascade)
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("IC MC y distribution")
	# ax.hist(energy_bins[:-1], energy_bins, weights = cascade / bin_widths,\
	# 				 label="ORCA cascade", histtype="step")
	# ax.hist(energy_bins[:-1], energy_bins, weights = track / bin_widths,\
	# 				 label="ORCA track", histtype="step")
	ax.hist(y_bins[:-1], y_bins, weights = ic_cascade / bin_widths,\
					 label="IC cascade", histtype="step")
	ax.hist(y_bins[:-1], y_bins, weights = ic_track / bin_widths,\
					 label="IC track", histtype="step")
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel("y")
	# ax.set_xlim(1, 50)
	ax.set_ylabel("event rate [3yrs]")
	# ax.grid(True)
	ax.legend()
	# plt.show()
	fig.savefig("./../ORCA/RecoPlots/IC_y_distribution")
	plt.close()

# y_distribution()

def track_E_vs_y_distribution():
	# ic_cas_weights = np.zeros_like(ic_sim.W_mc)
	ic_weights = np.zeros_like(ic_sim.W_mc)
	for i in range(len(ic_weights)):
		# ic_cas_weights[i] = ic_analysis.bf_weights[0][0][0][i] + ic_analysis.bf_weights[0][1][0][i] +  \
		# 					ic_analysis.bf_weights[1][0][0][i] + ic_analysis.bf_weights[1][1][0][i]
		ic_weights[i] = ic_analysis.bf_weights[0][0][1][i] + ic_analysis.bf_weights[0][1][1][i] +  \
							ic_analysis.bf_weights[1][0][1][i] + ic_analysis.bf_weights[1][1][1][i]	
	binnum = 10
	y_bins = np.logspace(-4, 0, binnum + 1)
	E_bins = np.logspace(0, 2, binnum + 1)
	y_bin_widths = np.zeros((binnum,))
	E_bin_widths = np.zeros((binnum,))
	for i in range(binnum):
		y_bin_widths[i] = y_bins[i+1] - y_bins[i]
		E_bin_widths[i] = E_bins[i+1] - E_bins[i]
	bin_areas = np.ndarray((binnum, binnum))
	for i in range(binnum):
		for j in range(binnum):
			bin_areas[i][j] = E_bin_widths[i] * y_bin_widths[j]
	# ic_cascade, _ = np.histogram(ic_sim.file["y"], bins = y_bins, weights = ic_cas_weights)
	ic_track, _, _ = np.histogram2d(ic_sim.E_tr, ic_sim.file["y"], bins = [E_bins, y_bins], weights = ic_weights)
	ic_track /= bin_areas
	# print(cascade)
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("ORCA MC energy vs y distribution")
	c = ax.pcolormesh(E_bins, y_bins, ic_track.T)
	plt.colorbar(c)
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel("Energy[GeV")
	ax.set_ylabel("y")
	ax.legend()
	# plt.show()
	fig.savefig("./../ORCA/RecoPlots/IC_track_E_y_distribution")
	plt.close()

# track_E_vs_y_distribution()

def cas_E_vs_y_distribution():
	# ic_cas_weights = np.zeros_like(ic_sim.W_mc)
	ic_weights = np.zeros_like(ic_sim.W_mc)
	for i in range(len(ic_weights)):
		ic_weights[i] = ic_analysis.bf_weights[0][0][0][i] + ic_analysis.bf_weights[0][1][0][i] +  \
							ic_analysis.bf_weights[1][0][0][i] + ic_analysis.bf_weights[1][1][0][i]
		# ic_weights[i] = ic_analysis.bf_weights[0][0][1][i] + ic_analysis.bf_weights[0][1][1][i] +  \
		# 					ic_analysis.bf_weights[1][0][1][i] + ic_analysis.bf_weights[1][1][1][i]	
	binnum = 10
	y_bins = np.logspace(-4, 0, binnum + 1)
	E_bins = np.logspace(0, 2, binnum + 1)
	y_bin_widths = np.zeros((binnum,))
	E_bin_widths = np.zeros((binnum,))
	for i in range(binnum):
		y_bin_widths[i] = y_bins[i+1] - y_bins[i]
		E_bin_widths[i] = E_bins[i+1] - E_bins[i]
	bin_areas = np.ndarray((binnum, binnum))
	for i in range(binnum):
		for j in range(binnum):
			bin_areas[i][j] = E_bin_widths[i] * y_bin_widths[j]
	# ic_cascade, _ = np.histogram(ic_sim.file["y"], bins = y_bins, weights = ic_cas_weights)
	ic_track, _, _ = np.histogram2d(ic_sim.E_tr, ic_sim.file["y"], bins = [E_bins, y_bins], weights = ic_weights)
	ic_track /= bin_areas
	# print(cascade)
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("IC MC cascades energy vs y distribution")
	c = ax.pcolormesh(E_bins, y_bins, ic_track.T)
	plt.colorbar(c)
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel("Energy [GeV]")
	ax.set_ylabel("y")
	ax.legend()
	# plt.show()
	fig.savefig("./../ORCA/RecoPlots/IC_cas_E__y_distribution")
	plt.close()

# cas_E_vs_y_distribution()