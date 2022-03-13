import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import util
from util import get_index
import scipy as scp
from scipy.interpolate import interp1d

class ICEffectiveAnalysis:
	def __init__(self, lo, hi, binnum, mc = 'neutrino_mc.csv', nusec = 'nu.txt', nubarsec = 'nubar.txt'):
		self.lo = lo 
		self.hi = hi 
		self.bins = np.logspace(np.log10(lo), np.log10(hi), num = binnum)
		self.widths = np.zeros(len(self.bins) - 1)
		for i in range(len(self.widths)):
			self.widths[i] = self.bins[i+1] - self.bins[i]
		self.nuxsection = pd.read_csv(nusec, sep = ' ', usecols = [0, 1, 2])
		self.nubarxsection = pd.read_csv(nubarsec, sep = ' ', usecols = [0, 1, 2])
		self.MC = pd.read_csv(mc)
		self.energy = self.MC["true_energy"]
		self.weights = self.MC["weight"]
		self.pdg = self.MC["pdg"]
		self.current = self.MC["current_type"]
		self.areas = np.zeros((2, 4, len(self.weights)))
		self.volumes = np.zeros((2, 4, len(self.weights)))

		self.nc_mask = self.MC["current_type"] == 0
		self.cc_mask = self.MC["current_type"] == 1
		self.nu_mask = (self.MC["pdg"]) / np.abs((self.MC["pdg"])) == 1
		self.nubar_mask = (self.MC["pdg"]) / np.abs((self.MC["pdg"])) == -1
		self.nue_mask = (self.MC["pdg"]) == 12
		self.numu_mask = (self.MC["pdg"]) == 14
		self.nutau_mask = (self.MC["pdg"]) == 16
		self.nuebar_mask = (self.MC["pdg"]) == -12
		self.numubar_mask = (self.MC["pdg"]) == -14
		self.nutaubar_mask = (self.MC["pdg"]) == -16

	def computeArea(self):
		for i in range(len(self.weights)):
			nutype, flavor = get_index(int(self.current[i]), int(self.pdg[i]))
			self.areas[nutype][flavor][i] = self.weights[i] / (4 * np.pi)

	def computeVolume(self):
		# first define the function that finds the cross sections
		findsec_nu = util.interpolate_xsection(1)
		findsec_nubar = util.interpolate_xsection(-1)
		# also define the number density
		nd = 0.9168 * (100 ** 3) * 6.022 * (10 ** 23)
		for i in range(len(self.weights)):
			# firt get the vector indices for later plugging in
			nutype, flavor = get_index(self.current[i], self.pdg[i])
			# now find the cross section of this event
			if self.energy[i] < 0.01:
				xsec = 0
			elif self.energy[i] > 125:
				xsec = 84.1363 + 26.4089
			else:
				if nutype == 0: # neutrino
					xsec = findsec_nu(self.energy[i]) * (10 ** -38) / (100 ** 2)
				elif nutype == 1:
					xsec = findsec_nubar(self.energy[i]) * (10 ** -38) / (100 ** 2)
			# now we can start filling in
			self.volumes[nutype][flavor][i] = self.weights[i] / (4* np.pi * xsec * nd)

	def plotArea(self):
		ehist, ebin = np.histogram(self.energy, bins = self.bins, weights = self.areas[0][0])
		ebarhist, ebarbin = np.histogram(self.energy, bins = self.bins, weights = self.areas[1][0])
		muhist, mubin = np.histogram(self.energy, bins = self.bins, weights = self.areas[0][1])
		mubarhist, mubarbin = np.histogram(self.energy, bins = self.bins, weights = self.areas[1][1])
		tauhist, taubin = np.histogram(self.energy, bins = self.bins, weights = self.areas[0][2])
		taubarhist, taubarbin = np.histogram(self.energy, bins = self.bins, weights = self.areas[1][2])
		nchist, ncbin = np.histogram(self.energy, bins = self.bins, weights = self.areas[0][3])
		ncbarhist, ncbarbin = np.histogram(self.energy, bins = self.bins, weights = self.areas[1][3])

		fig, ax = plt.subplots(figsize=(7,6))
		fig.suptitle("IceCube Effective Area")
		ax.hist(ebin[:-1], ebin, weights = ehist / self.widths, label=r"$\nu_eCC$", color="red", histtype="step")
		ax.hist(ebarbin[:-1], ebarbin, weights = ebarhist / self.widths, label=r"$\overline{\nu}_eCC$", \
								linestyle = '--', color="red", histtype="step")
		ax.hist(mubin[:-1], mubin, weights = muhist / self.widths, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
		ax.hist(mubarbin[:-1], mubarbin, weights = mubarhist / self.widths, label=r"$\overline{\nu}_{\mu}CC$", \
								linestyle = '--', color="blue", histtype="step")
		ax.hist(taubin[:-1], taubin, weights = tauhist / self.widths, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
		ax.hist(taubarbin[:-1], taubarbin, weights = taubarhist / self.widths, label=r"$\overline{\nu}_{\tau}CC$",\
								color="green", linestyle = "--", histtype="step")
		ax.hist(ncbin[:-1], ncbin, weights = nchist / self.widths, label=r"$\nu NC$", color="brown", histtype="step")
		ax.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist / self.widths, label=r"$\overline{\nu} NC$", color="brown",\
								linestyle = "--", histtype="step")

		ax.set_xscale("log")
		ax.set_yscale("log")
		ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
		ax.set_xlim(self.lo, self.hi)
		# ax.set_xlim(1, 50) # to check against ORCA 1-50 GeV plot
		ax.set_ylabel("Effective Area [m^2]")
		ax.grid(True)
		ax.legend(loc = 2)
		# plt.show()
		fig.savefig("./RecoPlots/IC_Effective_Area")
		plt.close()

	def plotVol(self):
		ehist, ebin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][0])
		ebarhist, ebarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][0])
		muhist, mubin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][1])
		mubarhist, mubarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][1])
		tauhist, taubin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][2])
		taubarhist, taubarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][2])
		nchist, ncbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][3])
		ncbarhist, ncbarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][3])

		fig, ax = plt.subplots(figsize=(7,6))
		fig.suptitle("IceCube Effective Volume")
		ax.hist(ebin[:-1], ebin, weights = ehist / self.widths, label=r"$\nu_eCC$", color="red", histtype="step")
		ax.hist(ebarbin[:-1], ebarbin, weights = ebarhist / self.widths, label=r"$\overline{\nu}_eCC$", \
								linestyle = '--', color="red", histtype="step")
		ax.hist(mubin[:-1], mubin, weights = muhist / self.widths, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
		ax.hist(mubarbin[:-1], mubarbin, weights = mubarhist / self.widths, label=r"$\overline{\nu}_{\mu}CC$", \
								linestyle = '--', color="blue", histtype="step")
		ax.hist(taubin[:-1], taubin, weights = tauhist / self.widths, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
		ax.hist(taubarbin[:-1], taubarbin, weights = taubarhist / self.widths, label=r"$\overline{\nu}_{\tau}CC$",\
								color="green", linestyle = "--", histtype="step")
		ax.hist(ncbin[:-1], ncbin, weights = nchist / self.widths, label=r"$\nu NC$", color="brown", histtype="step")
		ax.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist / self.widths, label=r"$\overline{\nu} NC$", color="brown",\
								linestyle = "--", histtype="step")
		ax.set_xscale("log")
		# ax.set_yscale("log")
		ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
		ax.set_xlim(self.lo, self.hi)
		ax.set_ylabel("Effective Volume [m^3]")
		ax.grid(True)
		ax.legend(loc = 2)
		# plt.show()
		fig.savefig("./RecoPlots/IC_Effective_volume")
		plt.close()

	def returnVol(self):
		ehist, ebin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][0])
		ebarhist, ebarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][0])
		muhist, mubin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][1])
		mubarhist, mubarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][1])
		tauhist, taubin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][2])
		taubarhist, taubarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][2])
		nchist, ncbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[0][3])
		ncbarhist, ncbarbin = np.histogram(self.energy, bins = self.bins, weights = self.volumes[1][3])
		ehist /= self.widths
		ebarhist /= self.widths
		muhist /= self.widths
		mubarhist /= self.widths
		tauhist /= self.widths
		taubarhist /= self.widths
		nchist /= self.widths
		ncbarhist /= self.widths
		return ehist, ebarhist, muhist, mubarhist, tauhist, taubarhist, nchist, ncbarhist
	
class ORCAEffectiveAnalysis:
	def __init__(self):
		# treat the histogram as "pseudoevents"
		self.evol = util.getORCAbins("./ORCA_Results/nueCC.csv")
		self.ebarvol = util.getORCAbins("./ORCA_Results/nuebarCC.csv")
		self.muvol = util.getORCAbins("./ORCA_Results/numuCC.csv")
		self.mubarvol = util.getORCAbins("./ORCA_Results/numubarCC.csv")
		self.tauvol = util.getORCAbins("./ORCA_Results/nutauCC.csv", tau = True)
		self.taubarvol = util.getORCAbins("./ORCA_Results/nutaubarCC.csv", tau = True)
		self.ncvol = util.getORCAbins("./ORCA_Results/nuNC.csv", nc = True)
		self.ncbarvol = util.getORCAbins("./ORCA_Results/nubarNC.csv", nc = True)
		self.earea = np.zeros(len(self.evol))
		self.ebararea = np.zeros(len(self.evol))
		self.muarea = np.zeros(len(self.evol))
		self.mubararea = np.zeros(len(self.evol))
		self.tauarea = np.zeros(len(self.evol))
		self.taubararea = np.zeros(len(self.evol))
		self.ncarea = np.zeros(len(self.evol))
		self.ncbararea = np.zeros(len(self.evol))
		self.e_bin = np.logspace(np.log10(1), np.log10(50), 51)
		self.e_center = np.zeros(len(self.e_bin) - 1)
		for i in range(len(self.e_center)):
			self.e_center[i] = (self.e_bin[i + 1] + self.e_bin[i]) / 2
		self.widths = np.zeros(len(self.e_center))
		for i in range(len(self.widths)):
			self.widths[i] = self.e_bin[i + 1] - self.e_bin[i]
		self.areas = np.zeros((2, 4, len(self.evol))) # This is not used for the moment

	# This multiplies each volume bin weight by the bin width, so the unit of volume is m^3 GeV
	def computeVolume(self):
		for i in range(len(self.evol)):
			self.evol[i] *= self.widths[i]
			self.muvol[i] *= self.widths[i]
			self.tauvol[i] *= self.widths[i]
			self.ncvol[i] *= self.widths[i]
			self.ebarvol[i] *= self.widths[i]
			self.mubarvol[i] *= self.widths[i]
			self.taubarvol[i] *= self.widths[i]
			self.ncbarvol[i] *= self.widths[i]
		
	# this multiplies the volume by xsection and number density, but not by 4pi, unit is m^2 GeV
	def computeArea(self):
		findsec_nu = util.interpolate_xsection(1)
		findsec_nubar = util.interpolate_xsection(-1)
		for i in range(len(self.evol)):
			# energy of the current pseudoevent
			energy = self.e_center[i]
			# find the cross section corresponding to nu and nubar respectively
			nusec = findsec_nu(energy) * (10 ** -38) / (100 ** 2)
			nubarsec = findsec_nubar(energy) * (10 ** -38) / (100 ** 2)
			# also define the number density
			nd = 0.9168 * (100 ** 3) * 6.022 * (10 ** 23)
			self.earea[i] = self.evol[i] * nd * nusec # * (4 * np.pi)
			self.muarea[i] = self.muvol[i] * nd * nusec # * (4 * np.pi)
			self.tauarea[i] = self.tauvol[i] * nd * nusec # * (4 * np.pi)
			self.ncarea[i] = self.ncvol[i] * nd * nusec # * (4 * np.pi)
			self.ebararea[i] = self.ebarvol[i] * nd * nubarsec # * (4 * np.pi)
			self.mubararea[i] = self.mubarvol[i] * nd * nubarsec # * (4 * np.pi)
			self.taubararea[i] = self.taubarvol[i] * nd * nubarsec # * (4 * np.pi)
			self.ncbararea[i] = self.ncbarvol[i] * nd * nubarsec # * (4 * np.pi)

	def plotVol(self):
		fig, ax = plt.subplots(figsize=(7,6))
		fig.suptitle("ORCA Effective Volume")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.evol / self.widths, label=r"$\nu_eCC$", color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ebarvol / self.widths, label=r"$\overline{\nu}_eCC$", \
								linestyle = '--', color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.muvol / self.widths, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.mubarvol / self.widths, label=r"$\overline{\nu}_{\mu}CC$", \
								linestyle = '--', color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.tauvol / self.widths, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.taubarvol / self.widths, label=r"$\overline{\nu}_{\tau}CC$",\
								color="green", linestyle = "--", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncvol / self.widths, label=r"$\nu NC$", color="brown", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncbarvol / self.widths, label=r"$\overline{\nu} NC$", color="brown",\
								linestyle = "--", histtype="step")
		ax.set_xscale("log")
		# ax.set_yscale("log")
		ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
		ax.set_ylabel("Effective Volume [m^3]")
		ax.set_xlim(1, 50)
		ax.grid(True)
		ax.legend(loc = 2)
		# plt.show()
		fig.savefig("./RecoPlots/ORCA_Effective_volume")
		plt.close()

	def plotArea(self):
		fig, ax = plt.subplots(figsize=(7,6))
		fig.suptitle("ORCA Effective Area")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.earea / self.widths, label=r"$\nu_eCC$", color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ebararea / self.widths, label=r"$\overline{\nu}_eCC$", \
								linestyle = '--', color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.muarea / self.widths, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.mubararea / self.widths, label=r"$\overline{\nu}_{\mu}CC$", \
								linestyle = '--', color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.tauarea / self.widths, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.taubararea / self.widths, label=r"$\overline{\nu}_{\tau}CC$",\
								color="green", linestyle = "--", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncarea / self.widths, label=r"$\nu NC$", color="brown", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncbararea / self.widths, label=r"$\overline{\nu} NC$", color="brown",\
								linestyle = "--", histtype="step")
		ax.set_xscale("log")
		ax.set_yscale("log")
		ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
		ax.set_ylabel("Effective Area [m^2]")
		ax.set_xlim(1, 50)
		ax.grid(True)
		ax.legend(loc = 2)
		# plt.show()
		fig.savefig("./RecoPlots/ORCA_Effective_Area")
		plt.close()

	def returnVol(self):
		evol = self.evol / self.widths
		ebarvol = self.ebarvol / self.widths
		muvol = self.muvol / self.widths
		mubarvol = self.mubarvol / self.widths
		tauvol = self.tauvol / self.widths
		taubarvol = self.taubarvol / self.widths
		ncvol = self.ncvol / self.widths
		ncbarvol = self.ncbarvol / self.widths

		return evol, ebarvol, muvol, mubarvol, tauvol, taubarvol, ncvol, ncbarvol



	# def make_csv(self):
	# 	#we make a neutrino_mc-like file that contains the orca information




# The followings are IC analysis
# ICanalysis = ICEffectiveAnalysis(1, 50, 51)
# ICanalysis.computeArea()
# ICanalysis.computeVolume()

# ORCAanalysis = ORCAEffectiveAnalysis()
# ORCAanalysis.computeArea()
# ORCAanalysis.computeVolume()

# ICe, ICeb, ICmu, ICmub, ICtau, ICtaub, ICnc, ICncb = ICanalysis.returnVol()
# ORe, OReb, ORmu, ORmub, ORtau, ORtaub, ORnc, ORncb = ORCAanalysis.returnVol()

# print(len(ICe))
# print(len(ORe))

# print(ICe / ORe)


