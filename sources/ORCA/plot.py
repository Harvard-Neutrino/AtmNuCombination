import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit
import math

import digitalizer as dgt
import util
from util import gaussian
from util import twogaussian
# import analyze as anl
from params import *

# input_file = pd.read_csv("ORCA.csv")
input_file = pd.read_csv(savename)
original = pd.read_csv("neutrino_mc.csv")
pid = input_file["pid"]
e_true = input_file["true_energy"]
e_reco = input_file["reco_energy"]
zen_true = input_file["true_zenith"]
zen_reco = input_file["reco_zenith"]

# for i in range(len(pid)):
#     if input_file["weight"][i] != 0:
#         if zen_reco[i] == -1:
#             print(i, " zen")
#         if e_reco[i] == -1:
#             print(i, " e")

def plot_mig_hist(binnum, top):

    tracks = dgt.Digitalizer(input_track, input_scale)
    cascades = dgt.Digitalizer(input_cascade, input_scale)

    tracks.set_palette(0, -3, 20)
    cascades.set_palette(0, -3, 20)

    tracks.digitalize(22, 22)
    cascades.digitalize(22, 22)

    D_tracks = tracks.extracted
    D_cascades = cascades.extracted

    if top == 1:
        data_entries = D_tracks[binnum]
    elif top == 0:
        data_entries = D_cascades[binnum]

    if not two_gaus:
        bins_centers = np.array([0.5 * (x_bins[i] + x_bins[i+1]) for i in range(len(x_bins)-1)])
        popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[binnum, 10, 1])

        # # Generate enough x values to make the curves look smooth.
        xspace = np.linspace(2, 55, 100000)

        # # Plot the histogram and the fitted function.
        plt.bar(bins_centers, data_entries, width=x_bins[1] - x_bins[0], color='navy', label=r'Histogram entries')
        plt.plot(xspace, gaussian(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')
    elif two_gaus:
        bins_centers = np.array([0.5 * (x_bins[i] + x_bins[i+1]) for i in range(len(x_bins)-1)])
        popt, pcov = curve_fit(twogaussian, xdata=bins_centers, ydata=data_entries, p0=[binnum, 10, 2, 1, .8])

        # # Generate enough x values to make the curves look smooth.
        xspace = np.linspace(2, 55, 100000)

        # # Plot the histogram and the fitted function.
        plt.bar(bins_centers, data_entries, width=x_bins[1] - x_bins[0], color='navy', label=r'Histogram entries')
        plt.plot(xspace, twogaussian(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')
    plt.xscale('log')
    if top == 1:
        plt.savefig("./MigMatPlots/Tracks/TrackMigMatBin{}.png".format(binnum))
    if top == 0:
        plt.savefig("./MigMatPlots/Cascades/CascadeMigMatBin{}.png".format(binnum))
    plt.close()

# for i in range(22):
#     plot_mig_hist(i, 0)
#     plot_mig_hist(i, 1)

# plots the zenith error line plot reproduction from ORCA paper
def plot_zenith_errors():
    e, ebar, mu, mubar = util.get_zenith_error()

    x = np.linspace(1.8, 53, 100000)
    ye = np.zeros_like(x)
    yeb = np.zeros_like(x)
    ymu = np.zeros_like(x)
    ymub = np.zeros_like(x)

    for i in range(len(x)):
        ye[i] = e(x[i])
        yeb[i] = ebar(x[i])
        ymu[i] = mu(x[i])
        ymub[i] = mubar(x[i])

    plt.plot(x, ye, color = "red", linewidth = 1.8, linestyle = '-', label = "Electron Neutrino")
    plt.plot(x, yeb, color = "red", linewidth = 1.8, linestyle = '--', label = "Electron AntiNeutrino")
    plt.plot(x, ymu, color = "blue", linewidth = 1.8, linestyle = '-', label = "Muon Neutrino")
    plt.plot(x, ymub, color = "blue", linewidth = 1.8, linestyle = '--', label = "Muon AntiNeutrino")

    plt.xscale('log')

    plt.xlabel("Neutrino Energy [Gev]")
    plt.ylabel("Median Angular Error [deg]")

    plt.title("KM3NeT Angular Error")
    plt.legend()
    plt.grid(True)

    plt.show()
    # plt.savefig("KM3NeT Angular Error")
    plt.close()

    return

# plot_zenith_errors()

# checks if digitalizer is working properly
def check_digitalizer_track():
    # first load the digitalized data
    tracks = dgt.Digitalizer(input_track, input_scale)
    cascades = dgt.Digitalizer(input_cascade, input_scale)

    tracks.set_palette(0, -3, 100)
    cascades.set_palette(0, -3, 100)

    tracks.digitalize(22, 22)
    cascades.digitalize(22, 22)

    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)

    Z = tracks.extracted
    Z[11][14] = 0.1
    for i in range(22):
        currcol = Z[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(1.85, 53)
    plt.ylim(1.85, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    # plt.savefig("./RecoPlots/ORCA_Reco_track_probability")
    plt.close()

# check_digitalizer_track()

# plot the fake ORCA MC energy reco track normalized prob
def plot_energy_reco_track_probability():
    track_mask = pid == 1
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true[track_mask], e_reco[track_mask], bins=(x, y), weights = input_file["weight"][track_mask])
    # attempt to manually normalize column
    for i in range(22):
        currcol = Z[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())

    plt.xlim(1.85, 53)
    plt.ylim(1.85, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.legend()
    plt.show()
    # plt.savefig("./RecoPlots/ORCA_Reco_track_probability_weighted")
    plt.close()

# plot_energy_reco_track_probability()

#########################################################NEEDS FURTHER IMPLEMENTATION########################################################
# plot the IC MC energy reco track normalized prob
def plot_IC_energy_reco_track_probability():
    IC_input_file = pd.read_csv("neutrino_mc.csv")
    pid = IC_input_file["pid"]
    e_true = IC_input_file["true_energy"]
    e_reco = IC_input_file["reco_energy"]
    zen_true = IC_input_file["true_zenith"]
    zen_reco = IC_input_file["reco_zenith"]
    track_mask = pid == 1
    fifteens = np.zeros((22,))
    eightyfives = np.zeros((22,))
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true[track_mask], e_reco[track_mask], bins=(x, y))
    # attempt to manually normalize column
    for i in range(22):
        currcol = Z[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
        
        # now also get the 15 and 85 percentiles
        emin = x[i]
        emax = x[i + 1]
        current_bin_recos = []
        for k in range(len(IC_input_file["true_energy"])):
            if IC_input_file["true_energy"][i] > emin and IC_input_file["true_energy"][i] <= emax:
                current_bin_recos.append(IC_input_file["true_energy"][i])
        length_of_current = len(current_bin_recos)
        sorted_current = length_of_current.sort()
        fifteens[i] = sorted_current[math.ceil(.15 * length_of_current)]
        eightyfives[i] = sorted_current[math.ceil(.85 * length_of_current)]


    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    plt.savefig("./RecoPlots/ICMC_Reco_track_probability")
    plt.close()

# plot_IC_energy_reco_track_probability()
###################################################################################################################

def plot_energy_reco_cascade_probability():
    cas_mask = pid == 0
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true[cas_mask], e_reco[cas_mask], bins=(x, y), weights = input_file["weight"][cas_mask])
    # attempt to manually normalize column
    for i in range(22):
        currcol = Z[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())

    plt.xlim(1.85, 53)
    plt.ylim(1.85, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.legend()
    plt.show()
    # plt.savefig("./RecoPlots/ORCA_Reco_cascade_probability_weighted")
    plt.close()


# plot_energy_reco_cascade_probability()

def plot_IC_energy_reco_cascade_probability():
    IC_input_file = pd.read_csv("neutrino_mc.csv")
    pid = IC_input_file["pid"]
    e_true = IC_input_file["true_energy"]
    e_reco = IC_input_file["reco_energy"]
    zen_true = IC_input_file["true_zenith"]
    zen_reco = IC_input_file["reco_zenith"]
    track_mask = pid == 1
    cas_mask = pid == 0
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true[cas_mask], e_reco[cas_mask], bins=(x, y))
    # attempt to manually normalize column
    for i in range(22):
        currcol = Z[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    plt.savefig("./RecoPlots/ICMC_Reco_cascade_probability")
    plt.close()

# plot_IC_energy_reco_cascade_probability()

def plot_energy_reco_intermediate_probability():
    int_mask = pid == 2
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true[int_mask], e_reco[int_mask], bins=(x, y))
    # attempt to manually normalize column
    for i in range(22):
        currcol = Z[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    plt.savefig("./RecoPlots/ORCA_Reco_intermediate_probability")
    plt.close()

# plot_energy_reco_intermediate_probability()

def plot_zenith_reco():
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(np.cos(zen_true), np.cos(zen_reco), bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.savefig("./RecoPlots/ORCA_coszen_Reco_new")
    plt.close()

# plot_zenith_reco()

def plot_IC_zenith_reco():
    IC_input_file = pd.read_csv("neutrino_mc.csv")
    pid = IC_input_file["pid"]
    e_true = IC_input_file["true_energy"]
    e_reco = IC_input_file["reco_energy"]
    zen_true = IC_input_file["true_zenith"]
    zen_reco = IC_input_file["reco_zenith"]
    track_mask = pid == 1
    cas_mask = pid == 0
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(np.cos(zen_true), np.cos(zen_reco), bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.savefig("./RecoPlots/IC_coszen_Reco")
    plt.close()

# plot_IC_zenith_reco()

def plot_zenith_reco_range(elo, ehi):
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    select_true = []
    select_reco = []
    for i in range(len(zen_true)):
        if e_true[i] >= elo and e_true[i] <= ehi:
            select_true.append(zen_true[i])
            select_reco.append(zen_reco[i])
    Z, xedges, yedges = np.histogram2d(np.cos(select_true), np.cos(select_reco), bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.savefig("./RecoPlots/ICMC_coszen_with_ORCA_Reco_in_range_{}_to_{}".format(elo, ehi))
    plt.close()

# plot_zenith_reco_range(1, 5)
# plot_zenith_reco_range(5, 15)
# plot_zenith_reco_range(15, 55)

# currently only plots the numu distributions
def plot_numu_CC_topology_distribution():
    # first define the bins
    bins = np.logspace(np.log10(2), np.log10(50), 31)
    one = np.ones((30,))
    numu_mask = input_file["pdg"] == 14
    cc_mask = input_file["current_type"] == 1
    cascade_mask = input_file["pid"] == 0
    track_mask = input_file["pid"] == 1
    intermediate_mask = input_file["pid"] == 2

    numu_CC_track_hist, _ = np.histogram(e_true[numu_mask & track_mask & cc_mask], bins = bins)
    numu_CC_cascade_hist, _ = np.histogram(e_true[numu_mask & cascade_mask & cc_mask], bins = bins)
    numu_CC_inter_hist, _ = np.histogram(e_true[numu_mask & intermediate_mask & cc_mask], bins = bins)
    numu_CC_all_hist, _ = np.histogram(e_true[numu_mask & cc_mask], bins = bins)

    print(numu_CC_track_hist)
    print(numu_CC_cascade_hist)
    print(numu_CC_inter_hist)
    print(numu_CC_all_hist)

    # exit(0)

    fig, ax = plt.subplots(figsize=(7,6))
    fig.suptitle("ORCA MC nu_mu CC topology fractions")
    # first the neutrino ones
    ax.hist(bins[:-1], bins, weights = numu_CC_track_hist / numu_CC_all_hist,\
                     label=r"$\nu_{\mu}$ CC Track", color="blue", histtype="step")
    ax.hist(bins[:-1], bins, weights = one - numu_CC_cascade_hist / numu_CC_all_hist,\
                     label=r"$\nu_{\mu}$ CC Track", color="red", histtype="step")

    ax.set_xscale("log")
    ax.set_xlabel("neutrino energy [GeV]")
    ax.set_xlim(2, 50)
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction")
    ax.grid(True)
    ax.legend(loc = 4)
    plt.show()
    # fig.savefig("./RecoPlots/unosc_numu_CC_Topology_Fraction")
    plt.close()

# plot_numu_CC_topology_distribution()
