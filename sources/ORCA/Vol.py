import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import util
import scipy as scp
from scipy.interpolate import interp1d

plot = True

def IC():
	hi = 125
	lo = 1.
	binnum = 61 # for 1-50 I use 51, for 1-125 I use 61 or 71

	# reads the xsection txt file
	nuxsection = pd.read_csv("nu.txt", sep = ' ', usecols = [0, 1, 2])
	nubarxsection = pd.read_csv("nubar.txt", sep = ' ', usecols = [0, 1, 2])

	# reads the ICMC file
	ICMC = pd.read_csv("neutrino_mc.csv")

	findsec_nu = util.interpolate_xsection(1)
	findsec_nubar = util.interpolate_xsection(-1)

	def findx(energy, nutype):
		# print("in findx")
		if energy < 0.01:
			return 2.538e-05 # PUT TO 0
		elif energy > 125:
			return 84.1363 + 26.4089
		else:
			if nutype == 1:
				return findsec_nu(energy) * (10 ** -38) / (100 ** 2)
			if nutype == -1:
				return findsec_nubar(energy) * (10 ** -38) / (100 ** 2)


	# initialize V_eff
	V_eff = np.zeros_like(ICMC["weight"])

	totV = 0
	eV = 0
	ebarV = 0
	muV = 0
	mubarV = 0
	tauV = 0
	taubarV = 0
	ncV = 0
	ncbarV = 0

	# initializa bin normalization
	energy_bins_fine = np.logspace(np.log10(lo), np.log10(hi), num=binnum)


	# now get the Veff function
	for i in range(len(ICMC["true_energy"])):
		# first calculate xsection
		x = findx(ICMC["true_energy"][i], ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]))
		# number of ice nucleon
		nd = 0.9168 * (100 ** 3) * 6.022 * (10 ** 23) #/ 18
		# translate from cm^3 to m^3
		V_eff[i] = (ICMC["weight"][i] / (x * nd))
		if ICMC["true_energy"][i] >= lo and ICMC["true_energy"][i] <= hi:
			totV += V_eff[i]
			if ICMC["current_type"][i] == 1:
				if np.abs(ICMC["pdg"][i]) == 12:
					if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
						eV += V_eff[i]
					else:
						ebarV += V_eff[i]
				elif np.abs(ICMC["pdg"][i]) == 14:
					if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
						muV += V_eff[i]
					else:
						mubarV += V_eff[i]
				elif np.abs(ICMC["pdg"][i]) == 16:
					if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
						tauV += V_eff[i]
					else:
						taubarV += V_eff[i]
			elif ICMC["current_type"][i] == 0:
				if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
					ncV += V_eff[i]
				else:
					ncbarV += V_eff[i]


	# print("total effective volume: ", totV)
	# print("e CC effective volume: ", eV)
	# print("ebar CC effective volume: ", ebarV)
	# print("mu CC effective volume: ", muV)
	# print("mubar CC effective volume: ", mubarV)
	# print("tau CC effective volume: ", tauV)
	# print("taubar CC effective volume: ", taubarV)
	# print("NC effective volume: ", ncV)
	# print("NC bar effective volume: ", ncbarV)

	nc_mask = ICMC["current_type"] == 0
	cc_mask = ICMC["current_type"] == 1
	nu_mask = (ICMC["pdg"]) / np.abs((ICMC["pdg"])) == 1
	nubar_mask = (ICMC["pdg"]) / np.abs((ICMC["pdg"])) == -1
	nue_mask = (ICMC["pdg"]) == 12
	numu_mask = (ICMC["pdg"]) == 14
	nutau_mask = (ICMC["pdg"]) == 16
	nuebar_mask = (ICMC["pdg"]) == -12
	numubar_mask = (ICMC["pdg"]) == -14
	nutaubar_mask = (ICMC["pdg"]) == -16

	print("start binning")
	# use numpy to bin the data
	ehist, ebin = np.histogram(ICMC["true_energy"][nue_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nue_mask & cc_mask], density = True)
	ebarhist, ebarbin = np.histogram(ICMC["true_energy"][nuebar_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nuebar_mask & cc_mask], density = True)
	muhist, mubin = np.histogram(ICMC["true_energy"][numu_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[numu_mask & cc_mask], density = True)
	mubarhist, mubarbin = np.histogram(ICMC["true_energy"][numubar_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[numubar_mask & cc_mask], density = True)
	tauhist, taubin = np.histogram(ICMC["true_energy"][nutau_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nutau_mask & cc_mask], density = True)
	taubarhist, taubarbin = np.histogram(ICMC["true_energy"][nutaubar_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nutaubar_mask & cc_mask], density = True)
	nchist, ncbin = np.histogram(ICMC["true_energy"][nc_mask & nu_mask], bins=energy_bins_fine, \
		weights=V_eff[nc_mask & nu_mask] / 3, density = True)
	ncbarhist, ncbarbin = np.histogram(ICMC["true_energy"][nc_mask & nubar_mask], bins=energy_bins_fine, \
		weights=V_eff[nc_mask & nubar_mask] / 3, density = True)
	
	if not plot:
		return ehist * eV, ebarhist * ebarV, muhist * muV, mubarhist * mubarV, tauhist * tauV, taubarhist * taubarV, nchist * ncV, ncbarhist * ncbarV


	# Now plot these numpy histograms
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("IceCube Effective Volume")

	# # with scaling
	# ax.hist(ebin[:-1], ebin, weights = ehist * eV, label=r"$\nu_eCC$", color="red", histtype="step")
	# ax.hist(ebarbin[:-1], ebarbin, weights = ebarhist * ebarV, label=r"$\overline{\nu}_eCC$", \
	# 						linestyle = '--', color="red", histtype="step")
	# ax.hist(mubin[:-1], mubin, weights = muhist * muV, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
	# ax.hist(mubarbin[:-1], mubarbin, weights = mubarhist * mubarV, label=r"$\overline{\nu}_{\mu}CC$", \
	# 						linestyle = '--', color="blue", histtype="step")
	# ax.hist(taubin[:-1], taubin, weights = tauhist * tauV, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
	# ax.hist(taubarbin[:-1], taubarbin, weights = taubarhist * taubarV, label=r"$\overline{\nu}_{\tau}CC$",\
	# 						color="green", linestyle = "--", histtype="step")
	# ax.hist(ncbin[:-1], ncbin, weights = nchist * ncV, label=r"$\nu NC$", color="brown", histtype="step")
	# ax.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist * ncbarV, label=r"$\overline{\nu} NC$", color="brown",\
	# 						linestyle = "--", histtype="step")

	# without scaling
	ax.hist(ebin[:-1], ebin, weights = ehist, label=r"$\nu_eCC$", color="red", histtype="step")
	ax.hist(ebarbin[:-1], ebarbin, weights = ebarhist, label=r"$\overline{\nu}_eCC$", \
							linestyle = '--', color="red", histtype="step")
	ax.hist(mubin[:-1], mubin, weights = muhist, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
	ax.hist(mubarbin[:-1], mubarbin, weights = mubarhist, label=r"$\overline{\nu}_{\mu}CC$", \
							linestyle = '--', color="blue", histtype="step")
	ax.hist(taubin[:-1], taubin, weights = tauhist, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
	ax.hist(taubarbin[:-1], taubarbin, weights = taubarhist, label=r"$\overline{\nu}_{\tau}CC$",\
							color="green", linestyle = "--", histtype="step")
	ax.hist(ncbin[:-1], ncbin, weights = nchist, label=r"$\nu NC$", color="brown", histtype="step")
	ax.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist, label=r"$\overline{\nu} NC$", color="brown",\
							linestyle = "--", histtype="step")


	ax.set_xscale("log")
	# ax.set_yscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
	ax.set_xlim(lo, hi)
	ax.set_ylabel("Effective Volume [m^3]")
	ax.grid(True)
	ax.legend(loc = 2)
	# plt.show()
	fig.savefig("./checkvol/density_noscale_method_1-125_60bin")

# IC()

def IC_vol_check():
	# this used to produce some numbers using the binned histogram to check
	hi = 50
	lo = 1.
	binnum = 31 # for 1-50 I use 51, for 1-125 I use 61 or 71

	# reads the xsection txt file
	nuxsection = pd.read_csv("nu.txt", sep = ' ', usecols = [0, 1, 2])
	nubarxsection = pd.read_csv("nubar.txt", sep = ' ', usecols = [0, 1, 2])

	# reads the ICMC file
	ICMC = pd.read_csv("neutrino_mc.csv")

	findsec_nu = util.interpolate_xsection(1)
	findsec_nubar = util.interpolate_xsection(-1)

	def findx(energy, nutype):
		# print("in findx")
		if energy < 0.01:
			return 2.538e-05 # PUT TO 0
		elif energy > 125:
			return 84.1363 + 26.4089
		else:
			if nutype == 1:
				return findsec_nu(energy) * (10 ** -38) / (100 ** 2)
			if nutype == -1:
				return findsec_nubar(energy) * (10 ** -38) / (100 ** 2)


	# initialize V_eff
	V_eff = np.zeros_like(ICMC["weight"])

	totV = 0
	eV = 0
	ebarV = 0
	muV = 0
	mubarV = 0
	tauV = 0
	taubarV = 0
	ncV = 0
	ncbarV = 0

	# initializa bin normalization
	energy_bins_fine = np.logspace(np.log10(lo), np.log10(hi), num=binnum)


	# now get the Veff function
	for i in range(len(ICMC["true_energy"])):
		# first calculate xsection
		x = findx(ICMC["true_energy"][i], ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]))
		# number of ice nucleon
		nd = 0.9168 * (100 ** 3) * 6.022 * (10 ** 23) #/ 18
		# translate from cm^3 to m^3
		V_eff[i] = (ICMC["weight"][i] / (x * nd))
		if ICMC["true_energy"][i] >= lo and ICMC["true_energy"][i] <= hi:
			totV += V_eff[i]
			if ICMC["current_type"][i] == 1:
				if np.abs(ICMC["pdg"][i]) == 12:
					if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
						eV += V_eff[i]
					else:
						ebarV += V_eff[i]
				elif np.abs(ICMC["pdg"][i]) == 14:
					if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
						muV += V_eff[i]
					else:
						mubarV += V_eff[i]
				elif np.abs(ICMC["pdg"][i]) == 16:
					if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
						tauV += V_eff[i]
					else:
						taubarV += V_eff[i]
			elif ICMC["current_type"][i] == 0:
				if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
					ncV += V_eff[i]
				else:
					ncbarV += V_eff[i]


	# print("total effective volume: ", totV)
	print("e CC effective volume: ", eV)
	print("ebar CC effective volume: ", ebarV)
	# print("mu CC effective volume: ", muV)
	# print("mubar CC effective volume: ", mubarV)
	# print("tau CC effective volume: ", tauV)
	# print("taubar CC effective volume: ", taubarV)
	# print("NC effective volume: ", ncV)
	# print("NC bar effective volume: ", ncbarV)

	nc_mask = ICMC["current_type"] == 0
	cc_mask = ICMC["current_type"] == 1
	nu_mask = (ICMC["pdg"]) / np.abs((ICMC["pdg"])) == 1
	nubar_mask = (ICMC["pdg"]) / np.abs((ICMC["pdg"])) == -1
	nue_mask = (ICMC["pdg"]) == 12
	numu_mask = (ICMC["pdg"]) == 14
	nutau_mask = (ICMC["pdg"]) == 16
	nuebar_mask = (ICMC["pdg"]) == -12
	numubar_mask = (ICMC["pdg"]) == -14
	nutaubar_mask = (ICMC["pdg"]) == -16

	print("start binning")
	# use numpy to bin the data
	ehist, ebin = np.histogram(ICMC["true_energy"][nue_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nue_mask & cc_mask], density = True)
	ebarhist, ebarbin = np.histogram(ICMC["true_energy"][nuebar_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nuebar_mask & cc_mask], density = True)
	muhist, mubin = np.histogram(ICMC["true_energy"][numu_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[numu_mask & cc_mask], density = True)
	mubarhist, mubarbin = np.histogram(ICMC["true_energy"][numubar_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[numubar_mask & cc_mask], density = True)
	tauhist, taubin = np.histogram(ICMC["true_energy"][nutau_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nutau_mask & cc_mask], density = True)
	taubarhist, taubarbin = np.histogram(ICMC["true_energy"][nutaubar_mask & cc_mask], bins=energy_bins_fine, \
		weights=V_eff[nutaubar_mask & cc_mask], density = True)
	nchist, ncbin = np.histogram(ICMC["true_energy"][nc_mask & nu_mask], bins=energy_bins_fine, \
		weights=V_eff[nc_mask & nu_mask] / 3, density = True)
	ncbarhist, ncbarbin = np.histogram(ICMC["true_energy"][nc_mask & nubar_mask], bins=energy_bins_fine, \
		weights=V_eff[nc_mask & nubar_mask] / 3, density = True)

	# now I have the bins, I can do some checks
	# First check the area under the histogram

	binwidth = np.zeros(len(ebin) - 1)
	for i in range(len(binwidth)):
		binwidth[i] = ebin[i + 1] - ebin[i]

	e_density_integral = 0
	rescaled_e_density = 0
	for i in range(len(ehist)):
		e_density_integral += ehist[i] * binwidth[i]
		rescaled_e_density += ehist[i] * binwidth[i] * eV

	ebar_density_integral = 0
	rescaled_ebar_density = 0
	for i in range(len(ebarhist)):
		ebar_density_integral += ebarhist[i] * binwidth[i]
		rescaled_ebar_density += ebarhist[i] * binwidth[i] * ebarV

	print("nu_e histogram integral area is, ", e_density_integral)
	print("nu_ebar histogram integral area is, ", ebar_density_integral)
	print("rescaled nu_e histogram integral area is, ", rescaled_e_density)
	print("rescaled nu_ebar histogram integral area is, ", rescaled_ebar_density)

	
# IC_vol_check()

def IC_area_check():

	inputfile = ICMC = pd.read_csv("neutrino_mc.csv")


	nc_mask = ICMC["current_type"] == 0
	cc_mask = ICMC["current_type"] == 1
	nu_mask = (ICMC["pdg"]) / np.abs((ICMC["pdg"])) == 1
	nubar_mask = (ICMC["pdg"]) / np.abs((ICMC["pdg"])) == -1
	nue_mask = (ICMC["pdg"]) == 12
	numu_mask = (ICMC["pdg"]) == 14
	nutau_mask = (ICMC["pdg"]) == 16
	nuebar_mask = (ICMC["pdg"]) == -12
	numubar_mask = (ICMC["pdg"]) == -14
	nutaubar_mask = (ICMC["pdg"]) == -16


	xe = np.zeros_like(inputfile["weight"][cc_mask & nue_mask])
	# xeb = np.zeros_like(inputfile["weight"][cc_mask & nuebar_mask])

	ye = np.zeros_like(inputfile["weight"][cc_mask & nue_mask])
	# yeb = np.zeros_like(inputfile["weight"][cc_mask & nuebar_mask])

	# print(len(xe), len(ye))
	# exit(0)

	for i in range(len(xe)):
		xe[i] = (inputfile["true_energy"][cc_mask & nue_mask][i])
		
		ye[i] = (inputfile["weight"][cc_mask & nue_mask][i])
	
	# for i in range(len(xeb)):
	# 	xeb[i] = (inputfile["true_energy"][cc_mask & nuebar_mask][i])
		
	# 	yeb[i] = (inputfile["weight"][cc_mask & nuebar_mask][i])

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("IC nu_e CC Effective Area Directly From MC")
	resfe = interp1d(xe, ye)
	# resfebar = interp1d(xeb, yeb)

	newx = np.arange(1, 100, 0.1)
	newye = resfe(newx)
	# newyebar = resfebar(newx)
	# ax.plot(newx, newye, '-', label = r"$\nu_eCC$")
	ax.plot(xe, ye)
	# ax.plot(newx, newyebar, '-', label = r"$\bar{\nu}_eCC$")
	ax.set_yscale("log")
	ax.set_xscale('log')
	# ax.legend(loc = 2)
	# plt.show()
	fig.savefig("./checkvol/ICMC Direct Effective Area point")

IC_area_check()



def ORCA():
	ehist = util.getORCAbins("./ORCA_Results/nueCC.csv")
	ebarhist = util.getORCAbins("./ORCA_Results/nuebarCC.csv")
	muhist = util.getORCAbins("./ORCA_Results/numuCC.csv")
	mubarhist = util.getORCAbins("./ORCA_Results/numubarCC.csv")
	tauhist = util.getORCAbins("./ORCA_Results/nutauCC.csv", tau = True)
	taubarhist = util.getORCAbins("./ORCA_Results/nutaubarCC.csv", tau = True)
	nchist = util.getORCAbins("./ORCA_Results/nuNC.csv")
	ncbarhist = util.getORCAbins("./ORCA_Results/nubarNC.csv")

	if not plot:
		return ehist, ebarhist, muhist, mubarhist, tauhist, taubarhist, nchist, ncbarhist

	energy_bins_fine = np.logspace(np.log10(1), np.log10(50), num=51)

	ebin = ebarbin = mubin = mubarbin = taubin = taubarbin = ncbin = ncbarbin = energy_bins_fine

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("ORCA Effective Volume")

	ax.hist(ebin[:-1], ebin, weights = ehist, label=r"$\nu_eCC$", color="red", histtype="step")
	ax.hist(ebarbin[:-1], ebarbin, weights = ebarhist, label=r"$\overline{\nu}_eCC$", \
							linestyle = '--', color="red", histtype="step")
	ax.hist(mubin[:-1], mubin, weights = muhist, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
	ax.hist(mubarbin[:-1], mubarbin, weights = mubarhist, label=r"$\overline{\nu}_{\mu}CC$", \
							linestyle = '--', color="blue", histtype="step")
	ax.hist(taubin[:-1], taubin, weights = tauhist, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
	ax.hist(taubarbin[:-1], taubarbin, weights = taubarhist, label=r"$\overline{\nu}_{\tau}CC$",\
							color="green", linestyle = "--", histtype="step")
	ax.hist(ncbin[:-1], ncbin, weights = nchist, label=r"$\nu NC$", color="brown", histtype="step")
	ax.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist, label=r"$\overline{\nu} NC$", color="brown",\
							linestyle = "--", histtype="step")


	ax.set_xscale("log")
	# ax.set_yscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
	ax.set_xlim(1, 50)
	ax.set_ylabel("Effective Volume [m^3]")
	ax.grid(True)
	ax.legend(loc = 2)
	# plt.show()
	fig.savefig("ORCA Effective Volume Reproduce")


def ratio():
	ICe, ICeb, ICmu, ICmub, ICtau, ICtaub, ICnc, ICncb = IC()
	ORe, OReb, ORmu, ORmub, ORtau, ORtaub, ORnc, ORncb = ORCA()

	ye = ORe / ICe 
	yeb = OReb / ICeb
	ymu = ORmu / ICmu
	ymub = ORmub / ICmub
	ync = ORnc / ICnc
	yncb = ORncb / ICncb
	
	volbins = np.logspace(np.log10(1), np.log10(50), num=51)
	x = volbins[:-1]

		

	
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("ORCA/IC Effective Volume Ratio")

	ax.plot(x, ye, label=r"$\nu_eCC$", color="red")
	ax.plot(x, yeb, label=r"$\overline{\nu}_eCC$", \
							linestyle = '--', color="red")
	ax.plot(x, ymu, label=r"$\nu_{\mu}CC$", color="blue")
	ax.plot(x, ymub, label=r"$\overline{\nu}_{\mu}CC$", \
							linestyle = '--', color="blue")
	ax.plot(x, ync, label=r"$\nu_eCC$", color="brown")
	ax.plot(x, yncb, label=r"$\overline{\nu}_eCC$", \
							linestyle = '--', color="brown")


	ax.set_xscale("log")
	# ax.set_yscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
	ax.set_xlim(1, 50)
	ax.set_ylabel("Effective Volume Ratio [m^3/m^3]")
	ax.grid(True)
	ax.legend(loc = 2)
	# plt.show()
	fig.savefig("ORCA_to_IC Effective Volume Ratio")

# ratio()



