import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq
import nuflux

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})

# Define path to file (you may need to change this to match your system)
input_file = "neutrino_mc.csv"

# Load the file using pandas
input_data = pd.read_csv(input_file)


# Define masks to identify different neutrino flavors
nue_mask = (np.abs(input_data["pdg"]) == 12)
numu_mask = (np.abs(input_data["pdg"]) == 14)
nutau_mask = (np.abs(input_data["pdg"]) == 16)

# Define masks to identify interaction species
cascade_mask = (np.abs(input_data["pid"]) == 0)
track_mask = (np.abs(input_data["pid"]) == 1)

# Define masks to identify different flavor/interaction combinations.
nc_mask = input_data["current_type"] == 0
cc_mask = input_data["current_type"] == 1
nue_cc_mask = nue_mask & cc_mask
numu_cc_mask = numu_mask & cc_mask
nutau_cc_mask = nutau_mask & cc_mask
nue_nc_mask = nue_mask & nc_mask
numu_nc_mask = numu_mask & nc_mask
nutau_nc_mask = nutau_mask & nc_mask 
nue_cascade_mask = nue_mask & cascade_mask
nue_track_mask = nue_mask & track_mask
nutau_cascade_mask = nutau_mask & cascade_mask
nutau_track_mask = nutau_mask & track_mask
numu_cascade_mask = numu_mask & cascade_mask
numu_track_mask = numu_mask & track_mask

# Define some energy bins (used throughout this notebook)
energy_bins_fine = np.logspace(0., 2., num=21)
energy_bins_course = np.logspace(0., 2., num=11)

units = nsq.Const()

interactions = False

# Set propagation bins
E_min = 10.0*units.GeV
E_max = 1.0e3*units.PeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3

# theta23 numeric values to probe sensitivity
t23min = 0.5
t23max = 1.1
t23step = 0.1
t23l = np.arange(t23min, t23max + t23step, t23step)

# Set the chi squared bins limits
E_bin_min = 10 * units.GeV
E_bin_max = 1000 * units.GeV
E_n_bins = 10

# Set up chi squared bins
bins = np.zeros_like(t23l.shape[0], E_n_bins)

# Set up parameters as in nu-fit5.0
theta12 = 0.185778
theta13 = 0.047611
theta23 = 0.272222







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

nsq_atm.Set_MixingAngle(0, 1, 0.185778)
nsq_atm.Set_MixingAngle(0, 2, 0.047611)
nsq_atm.Set_MixingAngle(1, 2, theta23)
nsq_atm.Set_SquaremassDifference(1, 7.42e-5)
nsq_atm.Set_SquaremassDifference(2, 2.517e-3)


nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
nsq_atm.EvolveState()

lifetime = 365*24*60*60
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
# Note that converting to mHz for the rate, as this is a more suitable unit for the IceCube Upgrade 
fig, ax = plt.subplots(figsize=(7,6))
fig.suptitle("Reconstructed Energy Rated Weight")
ax.hist(input_data["reco_energy"], bins=energy_bins_fine, \
      weights=input_data["rate_weight"][nue_cc_mask], \
      label=r"$\nu_{e,CC}$", color="blue", histtype="step")




