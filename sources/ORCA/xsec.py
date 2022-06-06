import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def calc_CCNC_total_xsec():
	write = True

	# these are the CC and NC xsec for e, mu, ebar, mubar
	nuxsection = pd.read_csv("./ORCA_Results/nu.txt", sep = ' ', usecols = [0, 1, 2]) 
	nubarxsection = pd.read_csv("./ORCA_Results/nubar.txt", sep = ' ', usecols = [0, 1, 2])

	nutauxsec = pd.read_csv("./ORCA_Results/nutau_xs.csv", sep = ',', header = None,  usecols = [0, 1])
	nutaubarxsec = pd.read_csv("./ORCA_Results/taubar_xs.csv", sep = ',', header = None,  usecols = [0, 1])

	res_e = nuxsection["Energy"]
	nuCC = nuxsection["sigmaCC"]
	nuNC = nuxsection["sigmaNC"]
	nubarCC = nubarxsection["sigmaCC"]
	nubarNC = nubarxsection["sigmaNC"]

	nutau_e = nutauxsec[0]
	nutau_xs = nutauxsec[1]
	nutaubar_e = nutaubarxsec[0]
	nutaubar_xs = nutaubarxsec[1]

	length = len(nutau_e)
	extracted = np.zeros(length)
	energies = np.zeros(length)
	for i in range(length):
		extracted[i] = nutau_xs[i]
		energies[i] = nutau_e[i]

	f_tau_xs = interp1d(energies, extracted)

	length = len(nutaubar_e)
	extracted = np.zeros(length)
	energies = np.zeros(length)
	for i in range(length):
		extracted[i] = nutaubar_xs[i]
		energies[i] = nutaubar_e[i]

	f_taubar_xs = interp1d(energies, extracted)


	# these are the CC nutau and nutaubar
	tau_ratio = pd.read_csv("./ORCA_Results/tau_xs_ratio.csv", header = None, usecols = [0, 1])

	length = len(tau_ratio[0])
	extracted = np.zeros(length)
	energies = np.zeros(length)
	for i in range(length):
		extracted[i] = tau_ratio[1][i]
		energies[i] = tau_ratio[0][i]

	f_tau_ratio = interp1d(energies, extracted)

	# print(f_tau_ratio(0))

	res_nu = np.zeros_like(res_e)
	res_nubar = np.zeros_like(res_e)
	res_tau = np.zeros_like(res_e)
	res_taubar = np.zeros_like(res_e)
	res_nc = np.zeros_like(res_e)
	res_ncbar = np.zeros_like(res_e)

	for i in range(len(res_e)):
		res_nu[i] = nuCC[i]
		res_nc[i] = nuNC[i]
		res_nubar[i] = nubarCC[i] 
		res_ncbar[i] = nubarNC[i]
		try:
			# res_tau[i] = nuCC[i] * f_tau_ratio(res_e[i])
			# res_taubar[i] = nubarCC[i] * f_tau_ratio(res_e[i])
			res_tau[i] = f_tau_xs(res_e[i])
			res_taubar[i] = f_taubar_xs(res_e[i]) # needs changing
		except: # this is for 100-125GeV events
			# print(i)
			# res_tau[i] = nuCC[i] * 0.74
			# res_taubar[i] = nubarCC[i] * 0.74
			res_tau[i] = 0
			res_taubar[i] = 0

	# now save to a csv file first
	if write:
		df_array = np.ndarray((7, len(res_e)))
		df_array[0] = res_e
		df_array[1] = res_nu
		df_array[2] = res_nubar
		df_array[3] = res_tau
		df_array[4] = res_taubar
		df_array[5] = res_nc 
		df_array[6] = res_ncbar
		df = pd.DataFrame(df_array.T, columns = ["energy", "nu", "nubar", "tau", "taubar", "nc", "ncbar"])
		df.to_csv("all_xsecs.csv")


	f_nuxsec = interp1d(res_e, res_nu)
	f_tauxsec = interp1d(res_e, res_tau)
	f_ncxsec = interp1d(res_e, res_nc)
	f_nubarxsec = interp1d(res_e, res_nubar)
	f_taubarxsec = interp1d(res_e, res_taubar)
	f_ncbarxsec = interp1d(res_e, res_ncbar)

	xval = np.linspace(1, 100, 1000)
	nuxsec = f_nuxsec(xval)
	tauxsec = f_tauxsec(xval)
	ncxsec = f_ncxsec(xval)
	nubarxsec = f_nubarxsec(xval)
	taubarxsec = f_taubarxsec(xval)
	ncbarxsec = f_ncbarxsec(xval)
	tauratio = f_tau_ratio(xval)

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("Cross Sections")

	ax.plot(xval, nuxsec, color = 'b', label = "CC numu")
	ax.plot(xval, tauxsec, color = 'g', label = "CC nutau")
	ax.plot(xval, ncxsec, color = 'brown', label = "NC")
	ax.plot(xval, nubarxsec, '--', color = 'b', label = "CC numubar")
	ax.plot(xval, taubarxsec, '--', color = 'g', label = "CC nutaubar")
	ax.plot(xval, ncbarxsec, '--', color = 'brown', label = "NCbar")
	# ax.plot(xval, tauratio, label = "tau ratio")

	ax.legend()

	# ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlim(1, 50)
	plt.show()

	return f_nuxsec, f_tauxsec

calc_CCNC_total_xsec()

def plot_xsec_from_mc():
	unsorted_mc = pd.read_csv("neutrino_mc.csv")
	mc = unsorted_mc.sort_values(by = ["true_energy"], ascending = True)
	energy = mc["true_energy"]
	area = mc["weight"] * 100 ** 2
	xsec = mc["xsec"]
	interaction = mc["interaction_type"]
	current = mc["current_type"]
	pdg = mc["pdg"]
	cc_mask = current == 1
	nc_mask = current == 0
	qe_mask = interaction == 0
	res_mask = interaction == 1
	dis_mask = interaction == 2
	nue_mask = pdg == 12
	nuebar_mask = pdg == -12
	numu_mask = pdg == 14
	nutau_mask = pdg == 16

	fig, axes = plt.subplots(3, 4, figsize = (30, 10))
	axes[0][0].plot(energy[nc_mask & qe_mask & nue_mask], xsec[nc_mask & qe_mask & nue_mask],\
												 '.', linewidth = 1, label = "nu_e NCQE xsections")
	axes[0][0].plot(energy[cc_mask & qe_mask & nue_mask], xsec[cc_mask & qe_mask & nue_mask],\
												 '.', linewidth = 1, label = "nu_e CCQE xsections")

	axes[0][1].plot(energy[nc_mask & res_mask & nue_mask], xsec[nc_mask & res_mask & nue_mask],\
												 '.', linewidth = 1, label = "nu_e NCRes xsections")
	axes[0][1].plot(energy[cc_mask & res_mask & nue_mask], xsec[cc_mask & res_mask & nue_mask],\
												 '.', linewidth = 1, label = "nu_e CCRes xsections")

	axes[0][2].plot(energy[nc_mask & dis_mask & nue_mask], xsec[nc_mask & dis_mask & nue_mask],\
												 '.', linewidth = 1, label = "nu_e NCDIS xsections")
	axes[0][2].plot(energy[cc_mask & dis_mask & nue_mask], xsec[cc_mask & dis_mask & nue_mask],\
												 '.', linewidth = 1, label = "nu_e CCDIS xsections")

	axes[1][0].plot(energy[nc_mask & qe_mask & numu_mask], xsec[nc_mask & qe_mask & numu_mask],\
												 '.', linewidth = 1, label = "nu_mu NCQE xsections")
	axes[1][0].plot(energy[cc_mask & qe_mask & numu_mask], xsec[cc_mask & qe_mask & numu_mask],\
												 '.', linewidth = 1, label = "nu_mu CCQE xsections")

	axes[1][1].plot(energy[nc_mask & res_mask & numu_mask], xsec[nc_mask & res_mask & numu_mask],\
												 '.', linewidth = 1, label = "nu_mu NCRes xsections")
	axes[1][1].plot(energy[cc_mask & res_mask & numu_mask], xsec[cc_mask & res_mask & numu_mask],\
												 '.', linewidth = 1, label = "nu_mu CCRes xsections")

	axes[1][2].plot(energy[nc_mask & dis_mask & numu_mask], xsec[nc_mask & dis_mask & numu_mask],\
												 '.', linewidth = 1, label = "nu_mu NCDIS xsections")
	axes[1][2].plot(energy[cc_mask & dis_mask & numu_mask], xsec[cc_mask & dis_mask & numu_mask],\
												 '.', linewidth = 1, label = "nu_mu CCDIS xsections")

	axes[2][0].plot(energy[nc_mask & qe_mask & nutau_mask], xsec[nc_mask & qe_mask & nutau_mask],\
												 '.', linewidth = 1, label = "nu_tau NCQE xsections")
	axes[2][0].plot(energy[cc_mask & qe_mask & nutau_mask], xsec[cc_mask & qe_mask & nutau_mask],\
												 '.', linewidth = 1, label = "nu_tau CCQE xsections")

	axes[2][1].plot(energy[nc_mask & res_mask & nutau_mask], xsec[nc_mask & res_mask & nutau_mask],\
												 '.', linewidth = 1, label = "nu_tau NCRes xsections")
	axes[2][1].plot(energy[cc_mask & res_mask & nutau_mask], xsec[cc_mask & res_mask & nutau_mask],\
												 '.', linewidth = 1, label = "nu_tau CCRes xsections")

	axes[2][2].plot(energy[nc_mask & dis_mask & nutau_mask], xsec[nc_mask & dis_mask & nutau_mask],\
												 '.', linewidth = 1, label = "nu_tau NCDIS xsections")
	axes[2][2].plot(energy[cc_mask & dis_mask & nutau_mask], xsec[cc_mask & dis_mask & nutau_mask],\
												 '.', linewidth = 1, label = "nu_tau CCDIS xsections")


	axes[0][0].set_ylabel("xsection [1e-38 cm^2], nu_e")
	axes[1][0].set_ylabel("xsection [1e-38 cm^2], nu_mu")
	axes[2][0].set_ylabel("xsection [1e-38 cm^2], nu_tau")

	axes[2][0].set_xlabel("Energy [GeV]")
	axes[2][1].set_xlabel("Energy [GeV]")
	axes[2][2].set_xlabel("Energy [GeV]")
	axes[2][3].set_xlabel("Energy [GeV]")

	axes[0][0].legend()
	axes[0][0].set_xlim(0, 40)
	# axes[0][0].set_ylim(1, 1.5)
	axes[0][0].set_yscale('log')

	axes[0][1].legend()
	axes[0][1].set_xlim(0, 40)
	# axes[0][1].set_ylim(1, 1.5)
	axes[0][1].set_yscale('log')

	axes[0][2].legend()
	axes[0][2].set_xlim(0, 100)
	axes[0][2].set_ylim(1, 150)
	axes[0][2].set_yscale('log')

	axes[1][0].legend()
	axes[1][0].set_xlim(0, 40)
	# axes[1][0].set_ylim(1, 0.2)
	axes[1][0].set_yscale('log')

	axes[1][1].legend()
	axes[1][1].set_xlim(0, 40)
	# axes[1][1].set_ylim(1, 0.2)
	axes[1][1].set_yscale('log')

	axes[1][2].legend()
	axes[1][2].set_xlim(0, 100)
	axes[1][2].set_ylim(1, 150)
	axes[1][2].set_yscale('log')

	axes[2][0].legend()
	axes[2][0].set_xlim(0, 40)
	# axes[2][0].set_ylim(1, 0.8)
	axes[2][0].set_yscale('log')

	axes[2][1].legend()
	axes[2][1].set_xlim(0, 40)
	# axes[2][1].set_ylim(1, 0.8)
	axes[2][1].set_yscale('log')

	axes[2][2].legend()
	axes[2][2].set_xlim(0, 100)
	axes[2][2].set_ylim(1, 150)
	axes[2][2].set_yscale('log')



	f_nuxsec, f_tauxsec = calc_CCNC_total_xsec()
	xval = np.linspace(1, 100, 1000)
	nuxsec = f_nuxsec(xval) # * (100 ** -2)
	tauxsec = f_tauxsec(xval) # * (100 ** -2)
	axes[0][3].plot(xval, nuxsec, label = "nu_e / nu_mu CC+NC xsections")
	axes[1][3].plot(xval, nuxsec, label = "nu_e / nu_mu CC+NC xsections")
	axes[2][3].plot(xval, tauxsec, label = "nu_tau CC+NC xsections")
	axes[0][3].legend()
	axes[1][3].legend()
	axes[2][3].legend()
	axes[0][3].set_yscale('log')
	axes[1][3].set_yscale('log')
	axes[2][3].set_yscale('log')

	plt.show()

# plot_xsec_from_mc()

