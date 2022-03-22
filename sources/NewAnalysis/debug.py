import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time
from matplotlib import pyplot as plt

import util
from params import *
import classdef as cl 



sim = cl.Simulation(pd.read_csv(input_file))
bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
fluxes = bf_fluxes
analysis = cl.Analysis(sim, bf_fluxes, fluxes)
util.get_all_weights(analysis, cl.PointType.BestFit)
util.get_all_weights(analysis, cl.PointType.Physical)

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

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("ORCA MC nu_{} {} topology fractions".format(flavname, curname))
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
	fig.savefig("./../ORCA/RecoPlots/nu{}_{}_Topology_Fraction".format(flavname, curname))
	plt.close()

ORCA_topology_details(0, 1)
ORCA_topology_details(1, 1)
ORCA_topology_details(2, 1)
ORCA_topology_details(3, 0)