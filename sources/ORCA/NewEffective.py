import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import util
from util import get_index
import scipy as scp
from scipy.interpolate import interp1d
from scipy.stats import linregress
import params


class ICEffectiveAnalysis:
	def __init__(self, lo, hi, binnum, mc = 'neutrino_mc.csv'):
		self.lo = lo 
		self.hi = hi 
		self.bins = np.logspace(np.log10(lo), np.log10(hi), num = binnum)
		self.widths = np.zeros(len(self.bins) - 1)
		for i in range(len(self.widths)):
			self.widths[i] = self.bins[i+1] - self.bins[i]
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

	# this only works for energies 0-100 GeV, but good enough for ORCA purposes
	def computeVolume(self):
		# first define the function that finds the cross sections
		# findsec_nu = util.interpolate_xsection(1)
		# findsec_nubar = util.interpolate_xsection(-1)

		# new implementation with tau xsec
		nusec, nubarsec, tausec, taubarsec, ncsec, ncbarsec = util.all_xsec()

		# findsec_nutau = util.interpolate_xsection(1, tau = True)
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
				if self.current[i] == 0:
					if nutype == 0:
						xsec = ncsec(self.energy[i]) * (10 ** -38) / (100 ** 2) * 3
					elif nutype == 1:
						xsec = ncbarsec(self.energy[i]) * (10 ** -38) / (100 ** 2) * 3
				elif self.current[1] == 1:
					if nutype == 0:
						if np.abs(self.pdg[i]) != 16:
							xsec = nusec(self.energy[i]) * (10 ** -38) / (100 ** 2)
						elif np.abs(self.pdg[i] == 16):
							xsec = tausec(self.energy[i]) * (10 ** -38) / (100 ** 2)
					elif nutype == 1:
						if np.abs(self.pdg[i]) != 16:
							xsec = nubarsec(self.energy[i]) * (10 ** -38) / (100 ** 2)
						elif np.abs(self.pdg[i]) == 16:
							xsec = taubarsec(self.energy[i]) * (10 ** -38) / (100 ** 2)

			# now we can start filling in
			if xsec != 0:
				self.volumes[nutype][flavor][i] = self.weights[i] / (4 * np.pi * xsec * nd)
			else:
				self.volumes[nutype][flavor][i] = 0

	# this does the computing in the order wanted
	def compute(self):
		self.computeArea()
		self.computeVolume()

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
		# ax.set_ylim(0, 2 * 1e6)
		# ax.set_yscale('log')
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
	def __init__(self, num_bins = 50):
		# treat the histogram as "pseudoevents"
		self.num_bins = num_bins
		self.evol = util.getORCAbins("./ORCA_Results/nueCC.csv")
		self.ebarvol = util.getORCAbins("./ORCA_Results/nuebarCC.csv")
		self.muvol = util.getORCAbins("./ORCA_Results/numuCC.csv")
		self.mubarvol = util.getORCAbins("./ORCA_Results/numubarCC.csv")
		self.tauvol = util.getORCAbins("./ORCA_Results/nutauCC.csv", tau = True)
		self.taubarvol = util.getORCAbins("./ORCA_Results/nutaubarCC.csv", tau = True)
		self.ncvol = util.getORCAbins("./ORCA_Results/nuNC.csv", nc = True)
		self.ncbarvol = util.getORCAbins("./ORCA_Results/nubarNC.csv", nc = True)
		self.earea = np.zeros(num_bins)
		self.ebararea = np.zeros(num_bins)
		self.muarea = np.zeros(num_bins)
		self.mubararea = np.zeros(num_bins)
		self.tauarea = np.zeros(num_bins)
		self.taubararea = np.zeros(num_bins)
		self.ncarea = np.zeros(num_bins)
		self.ncbararea = np.zeros(num_bins)
		self.e_bin = np.logspace(np.log10(1), np.log10(50), num_bins + 1)
		self.e_center = np.zeros(len(self.e_bin) - 1)
		for i in range(len(self.e_center)):
			self.e_center[i] = (self.e_bin[i + 1] + self.e_bin[i]) / 2
		self.widths = np.zeros(len(self.e_center))
		for i in range(len(self.widths)):
			self.widths[i] = self.e_bin[i + 1] - self.e_bin[i]

		self.original_e_bin = np.logspace(np.log10(1), np.log10(50), 50 + 1)
		self.original_e_center = np.zeros(len(self.original_e_bin) - 1)
		for i in range(len(self.original_e_center)):
			self.original_e_center[i] = (self.original_e_bin[i + 1] + self.original_e_bin[i]) / 2
		self.original_widths = np.zeros(len(self.original_e_center))
		for i in range(len(self.original_widths)):
			self.original_widths[i] = self.original_e_bin[i + 1] - self.original_e_bin[i]


		self.areas = np.zeros((2, 4, num_bins)) # This is not used for the moment

	# This multiplies each volume bin weight by the bin width, so the unit of volume is m^3 GeV
	def computeVolume(self):
		if self.num_bins != 50:
			original_bin = np.logspace(np.log10(1), np.log10(50), self.num_bins + 1)
			new_evol = np.zeros(self.num_bins)
			new_ebarvol = np.zeros(self.num_bins)
			new_muvol = np.zeros(self.num_bins)
			new_mubarvol = np.zeros(self.num_bins)
			new_tauvol = np.zeros(self.num_bins)
			new_taubarvol = np.zeros(self.num_bins)
			new_ncvol = np.zeros(self.num_bins)
			new_ncbarvol = np.zeros(self.num_bins)
			for i in range(len(self.e_center)):
				# print(i)
				energy = self.e_center[i]
				# print(energy)
				j = 0
				while j < len(self.original_e_center):
					if energy < self.original_e_center[j]:
						break
					j += 1
				# print(j)
				# print(self.original_e_center[j])
				new_evol[i] = self.evol[j]
				new_ebarvol[i] = self.ebarvol[j]
				new_muvol[i] = self.muvol[j]
				new_mubarvol[i] = self.mubarvol[j]
				new_tauvol[i] = self.tauvol[j]
				new_taubarvol[i] = self.taubarvol[j]
				new_ncvol[i] = self.ncvol[j]
				new_ncbarvol[i] = self.ncbarvol[j]
			self.evol = new_evol
			self.ebarvol = new_ebarvol
			self.muvol = new_muvol
			self.mubarvol = new_mubarvol
			self.tauvol = new_tauvol
			self.taubarvol = new_taubarvol
			self.ncvol = new_ncvol
			self.ncbarvol = new_ncbarvol


			
	# this multiplies the volume by xsection and number density, but not by 4pi, unit is m^2 GeV
	def computeArea(self):
		nusec, nubarsec, tausec, taubarsec, ncsec, ncbarsec = util.all_xsec()
		for i in range(len(self.earea)):
			# energy of the current pseudoevent
			energy = self.e_center[i]
			# find the cross section corresponding to nu and nubar respectively
			ccxsec = nusec(energy) * (10 ** -38) / (100 ** 2)
			ccbarxsec = nubarsec(energy) * (10 ** -38) / (100 ** 2)
			ncxsec = ncsec(energy) * (10 ** -38) / (100 ** 2)
			ncbarxsec = ncbarsec(energy) * (10 ** -38) / (100 ** 2) 
			tauxsec = tausec(energy)  * (10 ** -38) / (100 ** 2) 
			taubarxsec = taubarsec(energy)  * (10 ** -38) / (100 ** 2) 

			# also define the number density
			nd = 0.9168 * (100 ** 3) * 6.022 * (10 ** 23)

			self.earea[i] = self.evol[i] * nd * ccxsec # * (4 * np.pi)
			self.muarea[i] = self.muvol[i] * nd * ccxsec # * (4 * np.pi)
			self.tauarea[i] = self.tauvol[i] * nd * tauxsec # * (4 * np.pi)
			self.ncarea[i] = self.ncvol[i] * nd * ncxsec # * (4 * np.pi)
			self.ebararea[i] = self.ebarvol[i] * nd * ccbarxsec # * (4 * np.pi)
			self.mubararea[i] = self.mubarvol[i] * nd * ccbarxsec # * (4 * np.pi)
			self.taubararea[i] = self.taubarvol[i] * nd * taubarxsec # * (4 * np.pi)
			self.ncbararea[i] = self.ncbarvol[i] * nd * ncbarxsec # * (4 * np.pi)

	# this does the computing in the order wanted
	def compute(self):
		self.computeVolume()
		self.computeArea()

	def plotVol(self):
		fig, ax = plt.subplots(figsize=(7,6))
		fig.suptitle("ORCA Effective Volume")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.evol , label=r"$\nu_eCC$", color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ebarvol , label=r"$\overline{\nu}_eCC$", \
								linestyle = '--', color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.muvol , label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.mubarvol , label=r"$\overline{\nu}_{\mu}CC$", \
								linestyle = '--', color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.tauvol , label=r"$\nu_{\tau}CC$", color="green", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.taubarvol , label=r"$\overline{\nu}_{\tau}CC$",\
								color="green", linestyle = "--", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncvol , label=r"$\nu NC$", color="brown", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncbarvol , label=r"$\overline{\nu} NC$", color="brown",\
								linestyle = "--", histtype="step")
		ax.set_xscale("log")
		# ax.set_yscale("log")
		ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
		ax.set_ylabel("Effective Volume [m^3]")
		ax.set_xlim(1, 50)
		ax.grid(True)
		ax.legend(loc = 2)
		# plt.show()
		fig.savefig("./RecoPlots/ORCA_Effective_volume_20bins")
		plt.close()

	def plotArea(self):
		fig, ax = plt.subplots(figsize=(7,6))
		fig.suptitle("ORCA Effective Area")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.earea , label=r"$\nu_eCC$", color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ebararea , label=r"$\overline{\nu}_eCC$", \
								linestyle = '--', color="red", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.muarea , label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.mubararea , label=r"$\overline{\nu}_{\mu}CC$", \
								linestyle = '--', color="blue", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.tauarea , label=r"$\nu_{\tau}CC$", color="green", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.taubararea , label=r"$\overline{\nu}_{\tau}CC$",\
								color="green", linestyle = "--", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncarea , label=r"$\nu NC$", color="brown", histtype="step")
		ax.hist(self.e_bin[:-1], self.e_bin, weights = self.ncbararea , label=r"$\overline{\nu} NC$", color="brown",\
								linestyle = "--", histtype="step")
		ax.set_xscale("log")
		ax.set_yscale("log")
		ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
		ax.set_ylabel("Effective Area [m^2]")
		ax.set_xlim(1, 50)
		ax.grid(True)
		ax.legend(loc = 2)
		plt.show()
		# fig.savefig("./RecoPlots/ORCA_Effective_Area")
		plt.close()

	def returnVol(self):
		evol = self.evol #/ self.widths
		ebarvol = self.ebarvol #/ self.widths
		muvol = self.muvol #/ self.widths
		mubarvol = self.mubarvol #/ self.widths
		tauvol = self.tauvol #/ self.widths
		taubarvol = self.taubarvol #/ self.widths
		ncvol = self.ncvol #/ self.widths
		ncbarvol = self.ncbarvol #/ self.widths

		return evol, ebarvol, muvol, mubarvol, tauvol, taubarvol, ncvol, ncbarvol


def get_ratios(binnum = 20):
	ICanalysis = ICEffectiveAnalysis(1, 50, binnum + 1)
	ICanalysis.computeArea()
	ICanalysis.computeVolume()
	ORCAanalysis = ORCAEffectiveAnalysis(num_bins = binnum)
	ORCAanalysis.computeVolume()
	ORCAanalysis.computeArea()
	

	ICe, ICeb, ICmu, ICmub, ICtau, ICtaub, ICnc, ICncb = ICanalysis.returnVol()
	ORe, OReb, ORmu, ORmub, ORtau, ORtaub, ORnc, ORncb = ORCAanalysis.returnVol()

	# print(ICe / 1e6)
	# print(ORe / 1e6)

	# print(ICe.shape)

	energy_bins = np.logspace(np.log10(1), np.log10(50), binnum + 1)
	bin_centers = np.zeros(len(energy_bins) - 1)
	for i in range(len(bin_centers)):
		bin_centers[i] = (energy_bins[i + 1] + energy_bins[i]) / 2
	e_ratio = ORe / ICe
	eb_ratio = OReb / ICeb
	mu_ratio = ORmu / ICmu
	mub_ratio = ORmub / ICmub
	tau_ratio = ORtau / ICtau
	taub_ratio = ORtaub / ICtaub
	tau_ratio[6] = tau_ratio[7]
	taub_ratio[6] = taub_ratio[8]
	taub_ratio[7] = taub_ratio[8]
	nc_ratio = ORnc / ICnc
	ncb_ratio = ORncb / ICncb


	# f_e = linregress(np.log10(bin_centers), e_ratio)
	f_e = np.polyfit(np.log10(bin_centers)[:-1], e_ratio[:-1], deg = params.fit_deg)
	# f_mu = linregress(np.log10(bin_centers), mu_ratio)
	f_mu = np.polyfit(np.log10(bin_centers)[:-1], mu_ratio[:-1], deg = params.fit_deg)

	f_tau = linregress(np.log10(bin_centers)[6:], tau_ratio[6:])
	f_nc = linregress(np.log10(bin_centers), nc_ratio)

	# f_eb = linregress(np.log10(bin_centers), eb_ratio)
	f_eb = np.polyfit(np.log10(bin_centers)[:-1], eb_ratio[:-1], deg = params.fit_deg)

	# f_mub = linregress(np.log10(bin_centers), mub_ratio)
	f_mub = np.polyfit(np.log10(bin_centers)[:-1], mub_ratio[:-1], deg = params.fit_deg)

	f_taub = linregress(np.log10(bin_centers)[6:], taub_ratio[6:])
	f_ncb = linregress(np.log10(bin_centers), ncb_ratio)

	# print(e_A, e_b, e_k)

	xspace = np.logspace(np.log10(1), np.log10(50), 1000)

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("Effective Volume Ratios")

	ax.plot(bin_centers, e_ratio, color = 'r', linestyle = "-", alpha = 0.4)
	ax.plot(bin_centers, eb_ratio, color = 'r', linestyle = "--", alpha = 0.4)
	ax.plot(bin_centers, mu_ratio, color = 'b', linestyle = "-", alpha = 0.4)
	ax.plot(bin_centers, mub_ratio, color = 'b', linestyle = "--", alpha = 0.4)
	ax.plot(bin_centers, tau_ratio,  color = 'g', linestyle = "-", alpha = 0.4)
	ax.plot(bin_centers, taub_ratio, color = 'g', linestyle = "--", alpha = 0.4)
	ax.plot(bin_centers, nc_ratio, color = 'brown', linestyle = "-", alpha = 0.4)
	ax.plot(bin_centers, ncb_ratio, color = 'brown', linestyle = "--", alpha = 0.4)

	# ax.plot(xspace, f_e.intercept + f_e.slope * np.log10(xspace), label = "nu_e", color = 'r', linestyle = "-", alpha = 0.8)
	ax.plot(xspace, np.polyval(f_e, np.log10(xspace)), label = "nu_e", color = 'r', linestyle = "-", alpha = 0.8)

	# ax.plot(xspace, f_eb.intercept + f_eb.slope * np.log10(xspace), label = "nu_ebar", color = 'r', linestyle = "--", alpha = 0.8)
	ax.plot(xspace, np.polyval(f_eb, np.log10(xspace)), label = "nu_ebar", color = 'r', linestyle = "--", alpha = 0.8)

	# ax.plot(xspace, f_mu.intercept + f_mu.slope * np.log10(xspace), label = "nu_mu", color = 'b', linestyle = "-", alpha = 0.8)
	ax.plot(xspace, np.polyval(f_mu, np.log10(xspace)), label = "nu_mu", color = 'b', linestyle = "-", alpha = 0.8)

	# ax.plot(xspace, f_mub.intercept + f_mub.slope * np.log10(xspace), label = "nu_mubar", color = 'b', linestyle = "--", alpha = 0.8)
	ax.plot(xspace, np.polyval(f_mub, np.log10(xspace)), label = "nu_mubar", color = 'b', linestyle = "--", alpha = 0.8)

	ax.plot(xspace, f_tau.intercept + f_tau.slope * np.log10(xspace), label = "nu_tau", color = 'g', linestyle = "-", alpha = 0.8)
	ax.plot(xspace, f_taub.intercept + f_taub.slope * np.log10(xspace), label = "nu_taubar", color = 'g', linestyle = "--", alpha = 0.8)
	ax.plot(xspace, f_nc.intercept + f_nc.slope * np.log10(xspace), label = "nu_nc",color = 'brown', linestyle = "-", alpha = 0.8)
	ax.plot(xspace, f_ncb.intercept + f_ncb.slope * np.log10(xspace), label = "nu_ncbar",color = 'brown', linestyle = "--", alpha = 0.8)

	ax.set_ylim(0, 5)
	ax.set_xscale('log')

	# print(f_e.slope, f_e.intercept)
	# print(f_mu.slope, f_mu.intercept)
	# print(f_tau.slope, f_tau.intercept)
	# print(f_nc.slope, f_nc.intercept)

	ax.legend()
	# plt.show()
	# fig.savefig("./RecoPlots/Effective_Volume_Ratio_poly4_fit")

	return f_e, f_mu, f_tau, f_nc, f_eb, f_mub, f_taub, f_ncb

get_ratios(binnum = 20)

# ICanalysis = ICEffectiveAnalysis(1, 50, 21)
# ICanalysis.compute()
# ICanalysis.plotVol()
# ICanalysis.plotArea()

# ORanalysis = ORCAEffectiveAnalysis(num_bins = 20)
# ORanalysis.compute()
# ORanalysis.plotArea()



