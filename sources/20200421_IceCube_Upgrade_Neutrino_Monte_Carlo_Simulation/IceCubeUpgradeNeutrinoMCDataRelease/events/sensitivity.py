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

# List of theta13 and theta23 numeric values to probe
t13l = [.2, .3, .4, .5, .6, .7, .8]

# Prepare comparison values for flavor composition
nuel_13 = []
numul_13 = []
nutaul_13 = []

units = nsq.Const()

interactions = False

E_min = 10.0*units.GeV
E_max = 1.0e3*units.PeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3


for t13 in t13l:
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

    nsq_atm.Set_MixingAngle(0, 2, t13)
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
    nue = np.nansum(input_data["rate_weight"][nue_mask])
    numu = np.nansum(input_data["rate_weight"][numu_mask])
    nutau = np.nansum(input_data["rate_weight"][nutau_mask])
    nuel_13.append(nue)
    numul_13.append(numu)
    nutaul_13.append(nutau)
    
print(nuel_13)
print(numul_13)
print(nutaul_13)

chisq = [0, 0, 0, 0, 0, 0, 0]

def chisq(nuel, numul, nutaul):
    for i in range(0,6):
        res = (nuel[3] ** 2 - nuel[i] ** 2) / nuel[i] ** 2 \
                    + (numul[3] ** 2 - numul[i] ** 2) / numul[i] ** 2 \
                    + (nutaul[3] ** 2 - nutaul[i] ** 2) / nutaul[i] ** 2
        chisq[i] = res
   
    



