import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit

import digitalizer as dgt
import util
from util import gaussian
# import analyze as anl
from params import *

input_file = pd.read_csv("ORCA.csv")
pid = input_file["pid"]
e_true = input_file["true_energy"]
e_reco = input_file["reco_energy"]
zen_true = input_file["true_zenith"]
zen_reco = input_file["reco_zenith"]


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

    # if binnum == 11: # only manual hard code part
    #     data_entries[14] = 0.1

    # if binnum == 10:
    #     print(data_entries)
    #     data_entries[10] -= 0.34
    #     print(data_entries)

    bins_centers = np.array([0.5 * (x_bins[i] + x_bins[i+1]) for i in range(len(x_bins)-1)])
    popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[binnum, 10, 1])

    # # Generate enough x values to make the curves look smooth.
    xspace = np.linspace(2, 55, 100000)

    # # Plot the histogram and the fitted function.
    plt.bar(bins_centers, data_entries, width=x_bins[1] - x_bins[0], color='navy', label=r'Histogram entries')
    plt.plot(xspace, gaussian(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')
    if top == 1:
        plt.savefig("./MigMatPlots/Tracks/TrackMigMatBin{}.png".format(binnum))
    if top == 0:
        plt.savefig("./MigMatPlots/Cascades/CascadeMigMatBin{}.png".format(binnum))
    plt.close()

for i in range(22):
    plot_mig_hist(i, 0)
    plot_mig_hist(i, 1)

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

    # plt.show()
    plt.savefig("KM3NeT Angular Error")
    plt.close()

    return

# plot_zenith_errors()

def plot_energy_reco():
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true, e_reco, bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    plt.savefig("./RecoPlots/ICMC_with_ORCA_Reco")
    plt.close()

plot_energy_reco()

def plot_energy_reco_track():
    track_mask = pid == 1
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true[track_mask], e_reco[track_mask], bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    plt.savefig("./RecoPlots/ICMC_with_ORCA_Reco_track")
    plt.close()

plot_energy_reco_track()

def plot_energy_reco_cascade():
    cas_mask = pid == 0
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(e_true[cas_mask], e_reco[cas_mask], bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    plt.savefig("./RecoPlots/ICMC_with_ORCA_Reco_cas")
    plt.close()

plot_energy_reco_cascade()

def plot_zenith_reco():
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(np.cos(zen_true), np.cos(zen_reco), bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.savefig("./RecoPlots/ICMC_coszen_with_ORCA_Reco")
    plt.close()

plot_zenith_reco()

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

plot_zenith_reco_range(1, 5)
plot_zenith_reco_range(5, 15)
plot_zenith_reco_range(15, 55)
