import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux

from params import *

def propagate(theta23in, m31in, top):
    nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

    AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
    flux = nuflux.makeFlux('IPhonda2014_spl_solmin')
    for ic,cth in enumerate(nsq_atm.GetCosthRange()):
        for ie,E in enumerate(nsq_atm.GetERange()):
            nu_energy = E/units.GeV
            nu_cos_zenith = cth
            AtmInitialFlux[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
            AtmInitialFlux[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
            AtmInitialFlux[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
            AtmInitialFlux[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
            AtmInitialFlux[ic][ie][0][2] = 0 # flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau
            AtmInitialFlux[ic][ie][1][2] = 0 # flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar

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

        if input_data["true_energy"][i]*units.GeV > E_min and input_data["true_energy"][i]*units.GeV < E_max:
            rate_weight[i] = input_data["weight"][i]*nsq_atm.EvalFlavor(neuflavor,
                                                                        np.cos(input_data["true_zenith"][i]),
                                                                        input_data["true_energy"][i]*\
                                                                        units.GeV,neutype)*lifetime*meter_to_cm_sq*3 #3 years flux
            
    input_data["rate_weight"] = rate_weight
    if top == 0:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"][cascade_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][cascade_mask])
    elif top == 1:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"][track_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][track_mask])
    else:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"], bins = energy_bins_fine, weights = input_data["rate_weight"])


    return rate_weight, energy_hist_truth, energy_bins_truth