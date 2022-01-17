import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time

import util
from params import *
import classdef as cl 


def Me():
	NT23 = 24
	SqT23_min = 0.305
	SqT23_max = 0.705
	SqT23 = np.linspace(SqT23_min, SqT23_max, NT23+1, endpoint = True)
	SqT23BF = SqT23[17]
	T23 = np.arcsin(np.sqrt(SqT23))
	T23BF = T23[16]


	NDM31 = 24
	DM31_max = 3.0e-3
	DM31_min = 2.0e-3
	DM31 = np.linspace(DM31_min, DM31_max, NDM31+1, endpoint = True)
	DM31BF = DM31[12]

	NDCP = 19
	DCP_min = 0
	DCP_max = 2*np.pi
	DCP = np.linspace(DCP_min, DCP_max, NDCP+1, endpoint = True)
	DCPBF = DCP[13]
	
	sim = cl.Simulation(pd.read_csv(input_file))
	bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, T23BF, DM31BF, DCPBF)
	fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, np.arcsin(np.sqrt(0.50)), DM31BF, DCPBF)
	analysis = cl.Analysis(sim, bf_fluxes, fluxes)
	util.get_all_weights(analysis, cl.PointType.BestFit)
	util.get_all_weights(analysis, cl.PointType.Physical)

	# analysis.pre_apply_systematics(Syst)
	analysis.binning(cl.PointType.BestFit)
	analysis.binning(cl.PointType.Physical)

	Sys_BF = cl.Systematics(N_bf, delta_bf, gamma_bf, eps_bf, hv_bf, hv_bf)
	print(analysis.get_chisq(Sys_BF))

	# # now that we have the histogram, let's try calculating the chi squared due to numu cascade
	# numu_cas_bf = analysis.bf_histogram[0][1][0]
	# numu_cas_phy = analysis.histogram[0][1][0]

	# print(numu_cas_bf)
	# print(numu_cas_phy)

	# print(numu_cas_bf.shape)

	# # after we make sure this is the same as Ivan's right now we manually calculate chi squared
	# chisq = 0
	# for i in range(numu_cas_bf.shape[0]):
	# 	for j in range(numu_cas_bf.shape[1]):
	# 		if numu_cas_bf[i][j] > 0:
	# 			plus = (numu_cas_phy[i][j] - numu_cas_bf[i][j]) ** 2 / numu_cas_bf[i][j]
	# 			print(i, j, plus, "\n")
	# 			chisq += plus
	# print(chisq)

Me()

def Ivan():
	Time = 3.0*365*24*60*60
	meter_to_cm_sq = 1e4
	Radian = np.pi / 180.

	# Define path to file (you may need to change this to match your system)
	input_file = "neutrino_mc.csv"

	# Load the file using pandas
	input_data = pd.read_csv(input_file)

	NT23 = 24
	SqT23_min = 0.305
	SqT23_max = 0.705
	SqT23 = np.linspace(SqT23_min, SqT23_max, NT23+1, endpoint = True)
	SqT23BF = SqT23[17]
	T23 = np.arcsin(np.sqrt(SqT23))
	T23BF = T23[16]


	NDM31 = 24
	DM31_max = 3.0e-3
	DM31_min = 2.0e-3
	DM31 = np.linspace(DM31_min, DM31_max, NDM31+1, endpoint = True)
	DM31BF = DM31[12]

	NDCP = 19
	DCP_min = 0
	DCP_max = 2*np.pi
	DCP = np.linspace(DCP_min, DCP_max, NDCP+1, endpoint = True)
	DCPBF = DCP[13]

	# try flux propagation without classdef
	nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

	AtmIninue = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))

	flux = nuflux.makeFlux('IPhonda2014_spl_solmin')

	for ic,cth in enumerate(nsq_atm.GetCosthRange()):
		for ie,E in enumerate(nsq_atm.GetERange()):
			nu_energy = E/units.GeV
			nu_cos_zenith = cth 
			AtmIninue[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue

	WEBF = np.zeros((2, len(input_data["true_energy"])))

	nsq_atm.Set_MixingAngle(0,1, 33.44 * Radian)
	nsq_atm.Set_MixingAngle(0,2, 8.57 * Radian)
	nsq_atm.Set_MixingAngle(1,2, T23BF)
	nsq_atm.Set_SquareMassDifference(1,7.42e-5)
	nsq_atm.Set_SquareMassDifference(2,DM31BF)
	nsq_atm.Set_CPPhase(0,2,DCPBF)

	nsq_atm.Set_rel_error(1.0e-4)
	nsq_atm.Set_abs_error(1.0e-4)

	nsq_atm.Set_initial_state(AtmIninue,nsq.Basis.flavor)
	nsq_atm.EvolveState()

	for i in range(len(input_data["weight"])):
	#for i in range(len(rate_weight[mask_ErecMin])):

		if input_data["pdg"][i] > 0 :
			neutype = 0
		else:
			neutype = 1
					
		if input_data["pid"][i] == 0 :
			neuint = 0
		elif input_data["pid"][i] == 1 :
			neuint = 1

		if np.abs(input_data["pdg"][i]) == 12:
			neuflavor = 0
		elif np.abs(input_data["pdg"][i]) == 14:
			neuflavor = 1
		elif np.abs(input_data["pdg"][i]) == 16:
			neuflavor = 2


		#if input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:
		#WEBF[neuint][i] += input_data[mask_ErecMin]["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
		if input_data["reco_energy"][i] > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:
			WEBF[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

	print(WEBF[:10])