import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQUIDSpy as nsq
import nuflux
import seaborn as sns

from params import *

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})

# Obtain the rated weight of each event
# event topology is cascade 0 or track 1
def get_rated_weight_truth(top = 2):
    nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

    print("get_rated_weight_truth: propagating nu")
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

    nsq_atm.Set_MixingAngle(0, 1, theta12)
    nsq_atm.Set_MixingAngle(0, 2, theta13)
    nsq_atm.Set_MixingAngle(1, 2, theta23)
    nsq_atm.Set_SquareMassDifference(1, m21)
    nsq_atm.Set_SquareMassDifference(2, m31)

    nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
    nsq_atm.Set_ProgressBar(False) # progress bar will be printed on terminal
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
                                                                    units.GeV,neutype)*lifetime*meter_to_cm_sq*5

    # print("truth debug: before hist")
    input_data["rate_weight"] = rate_weight
    if top == 0:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"][cascade_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][cascade_mask])
    elif top == 1:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"][track_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][track_mask])
    else:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"], bins = energy_bins_fine, weights = input_data["rate_weight"])
    # print("truth debug: after hist")

    print("get_rated_weight_truth: energy rates: ", energy_hist_truth)

    return rate_weight , energy_hist_truth, energy_bins_truth


# Obtain binned energy given theta23 and m31 values
def get_energy_bins(theta23in, m31in, top = 2):
    nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

    AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
    flux = nuflux.makeFlux('honda2006')

    print("get_energy_bins: propagating nu")
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

    nsq_atm.Set_MixingAngle(0, 1, theta12)
    nsq_atm.Set_MixingAngle(0, 2, theta13)
    nsq_atm.Set_MixingAngle(1, 2, theta23in)
    nsq_atm.Set_SquareMassDifference(1, m21)
    nsq_atm.Set_SquareMassDifference(2, m31in)

    nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
    nsq_atm.Set_ProgressBar(False) # progress bar will be printed on terminal
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
                                                                    units.GeV,neutype)*lifetime*meter_to_cm_sq*5
    input_data["rate_weight"] = rate_weight
    
    # print("get_energy_bins_debug: before hist")
    # Now first obtain  the energy binned event rate distributions 1-100GeV
    if top == 0:
        energy_hist, energy_bins = np.histogram(input_data["reco_energy"][cascade_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][cascade_mask])
    elif top == 1:
        energy_hist, energy_bins = np.histogram(input_data["reco_energy"][track_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][track_mask])
    else:
        energy_hist, energy_bins = np.histogram(input_data["reco_energy"], bins = energy_bins_fine, weights = input_data["rate_weight"])
    # print("get_energy_bins_debug: after hist")
    print("get_energy_bins: energy rates: ", energy_hist)
    
    return energy_hist


# Get chisq for the contour plot
def get_chisq(t23, m31, top = 0):
    # Get the energy bins for the given t23, m31 and truth
    energy_bins = get_energy_bins(t23, m31, top)
    rate_weight_truth, energy_hist_truth, energy_bins_truth = get_rated_weight_truth(top)
    chisq = 0
    for i in range(len(energy_bins)):
        chisqplus = (energy_bins[i] - energy_hist_truth[i]) ** 2 /  energy_hist_truth[i]
        chisq += chisqplus
    return chisq

# Get the t23 chi sq raw profile (not minimizing over m31, set automatically to truth)
def get_t23_chi_profile(m31 = m31, top = 0):
    profile = np.zeros(len(t23l.tolist())).tolist()
    print(t23l)
    print(profile)
    for i in range(len(t23l.tolist())):
        profile[i] = get_chisq(t23l[i], m31, top)
        print(profile[i])
    return profile

# Get the m31 chi sq raw profile (not minimizing over t23, set automatically to truth)
def get_m31_chi_profile(t23 = theta23, top = 0):
    profile = np.zeros(len(m31l.tolist())).tolist()
    for i in range(len(m31l.tolist())):
        profile[i] = get_chisq(t23, m31l[i], top)
    return profile




# t23outfile = "t12bins.npy"
# m31outfile = "m31bins.npy"
# np.save(t23outfile, energy_hist_theta23)
# np.save(m31outfile, energy_hist_m31)



