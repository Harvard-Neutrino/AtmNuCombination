import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

hi = 50
lo = 1.

# reads the xsection txt file
xsection = pd.read_csv("xsec.txt", sep = ' ', usecols = [0, 1, 2])

# reads the ICMC file
ICMC = pd.read_csv("neutrino_mc.csv")

# def findx(energy, interaction_type):
# 	if interaction_type == 1: # CC
# 		amp = 84.1363
# 	elif interaction_type == 0: #NC
# 		amp = 26.4089
# 	for i in range(len(xsection["Energy"])):
# 		if energy <= xsection["Energy"][i]:
# 			if interaction_type == 1: # CC
# 				amp = xsection["sigmaCC"][i]
# 				break
# 			elif interaction_type == 0: #NC
# 				amp = xsection["sigmaNC"][i]
# 				break
# 			else:
# 				print("invalid interaction type detected")
# 				exit(1)

# 	return amp * (10 ** -38) / (100 ** 2)

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
# totcnt = 0
totV = 0

# initializa bin normalization
energy_bins_fine = np.logspace(np.log10(lo), np.log10(hi), num=51)
energy_bins_length = np.zeros((len(energy_bins_fine) - 1, ))
for i in range(len(energy_bins_length)):
	energy_bins_length[i] = (energy_bins_fine[i+1] - energy_bins_fine[i]) / 50

def findnorm(energy):
	norm = energy_bins_length[-1]
	for i in range(len(energy_bins_length)):
		if energy <= energy_bins_fine[i+1]:
			norm = energy_bins_length[i]
	return norm


# now get the Veff function
for i in range(len(ICMC["true_energy"])):
	# first calculate xsection
	x = findx(ICMC["true_energy"][i], ICMC["current_type"][i])
	# number of ice nucleon
	nd = 0.9168 * (100 ** 3) * 6.022 * (10 ** 23) / 18
	# get the normalization factor
	norm = findnorm(ICMC["true_energy"][i])
	# translate from cm^3 to m^3
	V_eff[i] = (ICMC["weight"][i] / (x * nd))
	if ICMC["true_energy"][i] >= lo and ICMC["true_energy"][i] <= hi:
		if ICMC["pdg"][i] == 12:
			if ICMC["current_type"][i] == 1:
				if ICMC["pdg"][i] / np.abs(ICMC["pdg"][i]) == 1:
					totV += V_eff[i]
	V_eff[i] *= (50 / norm)
	

# print(totcnt)
print(totV)

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



#now plot the effective volume
print("start plotting")
fig, ax = plt.subplots(figsize=(7,6))
fig.suptitle("IceCube Effective Volume")
# nu e 
ax.hist(ICMC["true_energy"][nue_mask & cc_mask], bins=energy_bins_fine, \
	  weights=V_eff[nue_mask & cc_mask], \
	  label=r"$\nu_eCC$", color="red", histtype="step")
ax.hist(ICMC["true_energy"][nuebar_mask & cc_mask], bins=energy_bins_fine, \
	  weights=V_eff[nuebar_mask & cc_mask], \
	  label=r"$\overline{\nu}_eCC$", color="red", linestyle = "--", histtype="step")
# nu mu
ax.hist(ICMC["true_energy"][numu_mask & cc_mask], bins=energy_bins_fine, \
	  weights=V_eff[numu_mask & cc_mask], \
	  label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
ax.hist(ICMC["true_energy"][numubar_mask & cc_mask], bins=energy_bins_fine, \
	  weights=V_eff[numubar_mask & cc_mask], \
	  label=r"$\overline{\nu}_{\mu}CC$", color="blue", linestyle = "--", histtype="step")
# nu tau
ax.hist(ICMC["true_energy"][nutau_mask & cc_mask], bins=energy_bins_fine, \
	  weights=V_eff[nutau_mask & cc_mask], \
	  label=r"$\nu_{\tau}CC$", color="green", histtype="step")
ax.hist(ICMC["true_energy"][nutaubar_mask & cc_mask], bins=energy_bins_fine, \
	  weights=V_eff[nutaubar_mask & cc_mask], \
	  label=r"$\overline{\nu}_{\tau}CC$", color="green", linestyle = "--", histtype="step")
# NC
ax.hist(ICMC["true_energy"][nc_mask & nu_mask], bins=energy_bins_fine, \
	  weights=V_eff[nc_mask & nu_mask] / 3, \
	  label=r"$\nu NC$", color="brown", histtype="step")
ax.hist(ICMC["true_energy"][nc_mask & nubar_mask], bins=energy_bins_fine, \
	  weights=V_eff[nc_mask & nubar_mask] / 3, \
	  label=r"$\overline{\nu} NC$", color="brown", linestyle = "--", histtype="step")
ax.set_xscale("log")
# ax.set_yscale("log")
ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
ax.set_xlim(lo, hi)
ax.set_ylabel("Effective Volume [m^3 / {}]".format(totV))
ax.grid(True)
ax.legend(loc = 2)
# fig.show()
fig.savefig("Effective Volume IC({} to {} GeV) nondensity norm 2.png".format(lo, hi))
plt.close()
