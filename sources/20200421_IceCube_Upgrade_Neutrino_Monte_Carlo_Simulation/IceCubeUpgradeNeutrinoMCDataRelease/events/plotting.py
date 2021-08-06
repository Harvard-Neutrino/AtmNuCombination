import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq
import nuflux
import seaborn as sns

from sensitivity import *

def distribution_plots():
    print("plotting sanity check distributions")
    
    # Plot the energy distribution
    fig, ax = plt.subplots(figsize=(7,6))
    fig.suptitle("Reco Energy Rated Distribution")
    ax.hist(input_data["reco_energy"][cascade_mask], bins=E_bin_plot, \
          weights=input_data["rate_weight"][cascade_mask], \
          label=r"$\nu_{All, Cascade}$", color="blue", histtype="step")
    ax.hist(input_data["reco_energy"][track_mask], bins=E_bin_plot, \
          weights=input_data["rate_weight"][track_mask], \
          label=r"$\nu_{All, Track}$", color="red", histtype="step")
    ax.set_xscale("log")
    ax.set_xlabel(r"$E_{\nu,\rm{reco}}$ [GeV]")
    ax.set_xlim(1, 100)
    ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
                     useOffset=None, useLocale=None, useMathText=None)
    ax.set_ylabel("Rate [5Years]")
    ax.grid(True)
    ax.legend()
    fig.savefig("Rate_For_Sensitivity.png", bbox_inches='tight')


    # Plot the angle distribution
    fig, ax = plt.subplots(figsize=(7,6))
    fig.suptitle("Reco Zenith Rated Distribution")
    ax.hist(np.cos(input_data["reco_zenith"][cascade_mask]), bins=cos_bin_plot, \
          weights=input_data["rate_weight"][cascade_mask], \
          label=r"$\nu_{All, Cascade}$", color="blue", histtype="step")
    ax.hist(np.cos(input_data["reco_zenith"][track_mask]), bins=cos_bin_plot, \
          weights=input_data["rate_weight"][track_mask], \
          label=r"$\nu_{All, Track}$", color="red", histtype="step")
    ax.set_xlabel(r"$\cos{\theta, \rm{reco}}$")
    ax.set_xlim(-1, 1)
    ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
                     useOffset=None, useLocale=None, useMathText=None)
    ax.set_ylabel("Rate [5Years]")
    ax.grid(True)
    ax.legend()
    fig.savefig("Zenith_Rate_For_Sensitivity.png", bbox_inches='tight')

    # Plot the overall distribution
    counts, _, _ = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=[energy_bins_fine, cos_bin_plot], \
          weights=input_data["rate_weight"])

    fig, ax = plt.subplots(figsize=(7,6))
    fig.suptitle("Reco Energy and Zenith Rated Distribution")
    ax.pcolormesh(energy_bins_fine, cos_bin_plot, counts.T)
    ax.set_xscale('log')
    ax.set_xlabel(r"$E_{\nu,\rm{reco}}$ [GeV]")
    ax.set_ylabel(r"$\cos{\theta, \rm{reco}}$")
    ax.set_xlim(1, 100)
    ax.set_ylim(-1, 1)
    ax.legend()
    fig.savefig("2D_Rate_For_Sensitivity.png", bbox_inches='tight')