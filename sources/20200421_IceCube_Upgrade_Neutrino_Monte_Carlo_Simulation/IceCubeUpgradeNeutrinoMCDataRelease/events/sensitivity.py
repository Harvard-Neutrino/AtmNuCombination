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
def get_rated_weight_truth():
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
    nsq_atm.Set_SquareMassDifference(1, 7.42e-5)
    nsq_atm.Set_SquareMassDifference(2, 2.517e-3)

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

    # Set the input rate weight and get the energy-binned fluxes for the ground truth
    input_data["rate_weight"] = get_rated_weight()
    energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"], bins = energy_bins_fine, weights = input_data["rate_weight"])

    return rate_weight, energy_hist_truth, energy_bins_truth

rate_weight_truth, energy_hist_truth, energy_bins_truth = get_rated_weight_truth()

# Obtain binned energy given theta23 and m31 values
def get_energy_bins(theta23in, m31in):
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
    
    # Now first obtain  the energy binned event rate distributions 1-100GeV
    energy_hist, energy_bins = np.histogram(input_data["reco_energy"], bins = energy_bins_fine, weights = input_data["rate_weight"])
    
    return energy_hist

t23sensitivity = True

if t23sensitivity:
    # Probe chi squared around truth value of theta23
    energy_hist_theta23 = np.zeros((len(t23l.tolist()), len(energy_bins_fine.tolist()) - 1)).tolist()
    for i in range(len(t23l.tolist())):
        print(i)
        energy_bins = get_energy_bins(t23l[i], m31).tolist()
        for j in range(len(energy_bins_fine.tolist()) - 1):
            print(j)  
            energy = energy_bins[j]
            energy_hist_theta23[i][j] = energy

    # Calculate non-normalized chi squared
    chisq = np.zeros((len(t23l.tolist()),))
    for i in range(len(t23l.tolist())):
        for j in range(len(energy_bins_fine.tolist())-1):
            chisqplus = (energy_hist_theta23[i][j] - energy_hist_truth[j]) ** 2 /  energy_hist_truth[j]
            chisq[i] += chisqplus


# plot un-normalized chisq for NH, probing values of t23
def plot_t23_chi():
    x = np.sin(t23l) ** 2
    y = chisq
    fig2, ax2 = plt.subplots(figsize=(7,6))
    fig2.suptitle("Chi-Sq NH")
    ax2.set_xlabel(r"$\sin^2{\theta_{23}}$")
    ax2.set_ylabel(r"$\chi^2_{NH}$")
    ax2.set_yscale("log")
    ax2.plot(x, y, color ="green")
    ax2.grid(True)
    fig2.savefig("t23_chi_sq(non-normal).png", bbox_inches='tight')

plot_t23_chi()


m31sensitivity = True

if m31sensitivity:
    # Probe chi squared around truth value of m31
    energy_hist_m31 = np.zeros((len(m31l.tolist()), len(energy_bins_fine.tolist()) - 1)).tolist()
    # print("energy_hist_theta23 initialization", energy_hist_m31)
    for i in range(len(m31l.tolist())):
        print(i)
        energy_bins = get_energy_bins(theta23, m31l[i]).tolist()
        for j in range(len(energy_bins_fine.tolist()) - 1):
            print(j)  
            energy = energy_bins[j]
            energy_hist_m31[i][j] = energy

    # Calculate non-normalized chi squared
    chisq2 = np.zeros((len(m31l.tolist()),))
    for i in range(len(m31l.tolist())):
        for j in range(len(energy_bins_fine.tolist())-1):
            chisqplus = (energy_hist_m31[i][j] - energy_hist_truth[j]) ** 2 /  energy_hist_truth[j]
            chisq2[i] += chisqplus

# plot un-normalized chisq for NH, probing values of m31
def plot_m31_chi():
    x = m31l
    y = chisq2
    fig3, ax3 = plt.subplots(figsize=(7,6))
    fig3.suptitle("Chi-Sq NH")
    ax3.set_xlabel(r"$m^2_{31}$")
    ax3.set_ylabel(r"$\sin^2{\theta_{23}}$")
    ax3.plot(x, y, color ="green")
    ax3.set_yscale("log")
    ax3.grid(True)
    fig3.savefig("m31_chi_sq(non-normal).png", bbox_inches='tight')
plot_m31_chi()




# t23outfile = "t12bins.npy"
# m31outfile = "m31bins.npy"
# np.save(t23outfile, energy_hist_theta23)
# np.save(m31outfile, energy_hist_m31)

# Get chisq for the contour plot
def get_chisq(t23, m31):
    energy_bins = get_energy_bins(t23, m31)
    chisq = 0
    for i in len(energy_bins):
        chisqplus = (energy_bins[i] - energy_hist_truth[i]) ** 2 /  energy_hist_truth[i]
        chisq += chisqplus
    return chisq


def plot_contour_chi():
    t23 = np.arange(t23min, t23max + t23step, t23step)
    m31 = np.arange(m31min, m31max + m31step, m31step)

    X, Y = np.meshgrid(np.sin(t23) ** 2, m31)
    Z = get_chisq(X, Y)

    fig4, ax4  = plt.subplots(figsize=(7,6))
    fig4.suptitle("Chi-Sq Contour NH")
    ax3.set_xlabel(r"$\sin^2{\theta_{23}}$")
    ax4.set_ylabel(r"$m^2_{31}$")
    axim = ax4.contourf(X,Y,Z,levels=[1e0,1e1,1e2,1e3],cmap=plt.cm.jet,norm = LogNorm())
    cb   = fig.colorbar(axim)

    fig4.savefig("Chisq Contour.png", bbox_inches="tight")

plot_contour_chi()