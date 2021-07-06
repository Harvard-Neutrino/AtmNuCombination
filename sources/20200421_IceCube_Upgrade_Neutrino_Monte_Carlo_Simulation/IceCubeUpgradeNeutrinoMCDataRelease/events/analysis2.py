import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQuIDSpy as nsq
import nuflux

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})


# Define masks to identify different neutrino flavors
nue_mask = (np.abs(input_data["pdg"]) == 12)
numu_mask = (np.abs(input_data["pdg"]) == 14)
nutau_mask = (np.abs(input_data["pdg"]) == 16)

# Define masks to identify different flavor/interaction combinations.
nc_mask = input_data["current_type"] == 0
cc_mask = input_data["current_type"] == 1
nue_cc_mask = nue_mask & cc_mask
numu_cc_mask = numu_mask & cc_mask
nutau_cc_mask = nutau_mask & cc_mask 



# Define path to file (you may need to change this to match your system)
input_file = "neutrino_mc.csv"

# Load the file using pandas
input_data = pd.read_csv(input_file)

units = nsq.Const()

interactions = False

E_min = 1.0*units.GeV
E_max = 1.0e3*units.PeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3

nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
flux = nuflux.makeFlux('honda2006')
for ic,cth in enumerate(nsq_atm.GetCosthRange()):
    for ie,E in enumerate(nsq_atm.GetERange()):
        nu_energy = E/units.GeV
        nu_cos_zenith = cth
        AtmInitialFlux[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
        AtmInitialFlux[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
        AtmInitialFlux[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
        AtmInitialFlux[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
        AtmInitialFlux[ic][ie][0][2] = flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau
        AtmInitialFlux[ic][ie][1][2] = flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar


nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
nsq_atm.EvolveState()

lifetime = 5.0*365*24*60*60
meter_to_cm_sq = 1e4
rate_weight = np.zeros_like(input_data["weight"])
for i in range(len(rate_weight)):
    if input_data["pdg"][i] > 0 :
        neutype = 0
    else:
        neutype = 1
        
    if np.abs(input_data["pdg"][i]) == 12:
        neuflavor = 0
    elif np.abs(input_data["pdg"][i]) == 14:
        neuflavor = 1
    elif np.abs(input_data["pdg"][i]) == 16:
        neuflavor = 2
        
    if input_data["true_energy"][i]*units.GeV < E_min or input_data["true_energy"][i]*units.GeV > E_max:
        rate_weight[i] = 0
        continue
    rate_weight[i] = input_data["weight"][i]*nsq_atm.EvalFlavor(neuflavor,
                                                                np.cos(input_data["true_zenith"][i]),
                                                                input_data["true_energy"][i]*\
                                                                units.GeV,neutype)*lifetime*meter_to_cm_sq
input_data["rate_weight"] = rate_weight

def plot_rate():
    # Note that converting to mHz for the rate, as this is a more suitable unit for the IceCube Upgrade 
    fig, ax = plt.subplots(figsize=(7,6))
    ax.hist(input_data["true_energy"][nue_cc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][nue_cc_mask], \
          label=r"$\nu_{e,CC}$", color="blue", histtype="step")
    ax.hist(input_data["true_energy"][numu_cc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][numu_cc_mask], \
          label=r"$\nu_{\mu,CC}$", color="red", histtype="step")
    ax.hist(input_data["true_energy"][nutau_cc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][nutau_cc_mask], \
          label=r"$\nu_{\tau,CC}$", color="green", histtype="step")
    ax.hist(input_data["true_energy"][nc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][nc_mask], \
          label=r"$\nu_{NC}$", color="grey", histtype="step")
    ax.set_xscale("log")
    ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
    ax.set_ylabel("Rate [mHz]")
    ax.grid(True)
    _ = ax.legend()
    fig.savefig("rated_weight_distribution_true_energy_2.png")

plot_rate()












