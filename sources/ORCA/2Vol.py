import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

hi = 50
lo = 1.

# reads the xsection txt file
xsection = pd.read_csv("xsec.txt", sep = ' ', usecols = [0, 1, 2])

# reads the ICMC file
ICMC = pd.read_csv("neutrino_mc.csv")

def findx(energy, interaction_type):
	amp = 84.1363+26.4089
	for i in range(len(xsection["Energy"])):
		if energy <= xsection["Energy"][i]:
			if interaction_type == 1 or interaction_type == 0:
				amp = xsection["sigmaCC"][i] + xsection["sigmaNC"][i]
				break
			else:
				print("invalid interaction type detected")
				exit(1)

	return amp * (10 ** -38) / (100 ** 2)

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
energy_bins_fine = np.logspace(np.log10(lo), np.log10(hi), num=31)


# now get the Veff function
for i in range(len(ICMC["true_energy"])):
	# first calculate xsection
	x = findx(ICMC["true_energy"][i], ICMC["current_type"][i])
	# number of ice nucleon
	nd = 0.9168 * (100 ** 3) * 6.022 * (10 ** 23) / 18
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


print("total effective volume: ", totV)
print("e CC effective volume: ", eV)
print("ebar CC effective volume: ", ebarV)
print("mu CC effective volume: ", muV)
print("mubar CC effective volume: ", mubarV)
print("tau CC effective volume: ", tauV)
print("taubar CC effective volume: ", taubarV)
print("NC effective volume: ", ncV)
print("NC bar effective volume: ", ncbarV)

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


# Now plot these numpy histograms
fig, ax = plt.subplots(figsize=(7,6))
fig.suptitle("IceCube Effective Volume")

ax.hist(ebin[:-1], ebin, weights = ehist * eV, label=r"$\nu_eCC$", color="red", histtype="step")
ax.hist(ebarbin[:-1], ebarbin, weights = ebarhist * ebarV, label=r"$\overline{\nu}_eCC$", \
						linestyle = '--', color="red", histtype="step")
ax.hist(mubin[:-1], mubin, weights = muhist * muV, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
ax.hist(mubarbin[:-1], mubarbin, weights = mubarhist * mubarV, label=r"$\overline{\nu}_{\mu}CC$", \
						linestyle = '--', color="blue", histtype="step")
ax.hist(taubin[:-1], taubin, weights = tauhist * tauV, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
ax.hist(taubarbin[:-1], taubarbin, weights = taubarhist * taubarV, label=r"$\overline{\nu}_{\tau}CC$",\
						 color="green", linestyle = "--", histtype="step")
ax.hist(ncbin[:-1], ncbin, weights = nchist * ncV, label=r"$\nu NC$", color="brown", histtype="step")
ax.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist * ncbarV, label=r"$\overline{\nu} NC$", color="brown",\
						 linestyle = "--", histtype="step")


ax.set_xscale("log")
# ax.set_yscale("log")
ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
ax.set_xlim(lo, hi)
ax.set_ylabel("Effective Volume [m^3]")
ax.grid(True)
ax.legend(loc = 2)
plt.show()
# fig.savefig("IC Effective Volume Normalized")

