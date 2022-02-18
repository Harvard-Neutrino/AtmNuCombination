import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

lo = 1
hi = 1000

# reads the ICMC file
ICMC = pd.read_csv("neutrino_mc.csv")

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
energy_bins_fine = np.logspace(np.log10(lo), np.log10(hi), num=51)

totA = 0
for i in range(len(ICMC["true_energy"])):
	if ICMC["true_energy"][i] >= lo and ICMC["true_energy"][i] <= hi:
		totA += ICMC["weight"][i]

#now plot the effective volume
print("start plotting")
fig, ax = plt.subplots(figsize=(7,6))
fig.suptitle("IceCube Effective Area")
# nu e 
ax.hist(ICMC["true_energy"][nue_mask & cc_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][nue_mask & cc_mask], density = True, \
	  label=r"$\nu_eCC$", color="red", histtype="step")
ax.hist(ICMC["true_energy"][nuebar_mask & cc_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][nuebar_mask & cc_mask], density = True, \
	  label=r"$\overline{\nu}_eCC$", color="red", linestyle = "--", histtype="step")
# nu mu
ax.hist(ICMC["true_energy"][numu_mask & cc_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][numu_mask & cc_mask], density = True, \
	  label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
ax.hist(ICMC["true_energy"][numubar_mask & cc_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][numubar_mask & cc_mask], density = True, \
	  label=r"$\overline{\nu}_{\mu}CC$", color="blue", linestyle = "--", histtype="step")
# nu tau
ax.hist(ICMC["true_energy"][nutau_mask & cc_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][nutau_mask & cc_mask], density = True, \
	  label=r"$\nu_{\tau}CC$", color="green", histtype="step")
ax.hist(ICMC["true_energy"][nutaubar_mask & cc_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][nutaubar_mask & cc_mask], density = True, \
	  label=r"$\overline{\nu}_{\tau}CC$", color="green", linestyle = "--", histtype="step")
# NC
ax.hist(ICMC["true_energy"][nc_mask & nu_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][nc_mask & nu_mask] / 3, density = True, \
	  label=r"$\nu NC$", color="brown", histtype="step")
ax.hist(ICMC["true_energy"][nc_mask & nubar_mask], bins=energy_bins_fine, \
	  weights=ICMC["weight"][nc_mask & nubar_mask] / 3, density = True, \
	  label=r"$\overline{\nu} NC$", color="brown", linestyle = "--", histtype="step")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
ax.set_xlim(lo, hi)
ax.set_ylabel("Effective Area [m^2] / {}".format(totA))
ax.grid(True)
ax.legend(loc = 2)
# fig.show()
fig.savefig("Effective Area IC ({}-{}GeV)".format(lo, hi))
plt.close()
