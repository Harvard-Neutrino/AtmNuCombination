import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit

import digitalizer as dgt
import util
from util import gaussian
import analyze as anl
from params import *

def plot_mig_hist(binnum):

    tracks = dgt.Digitalizer(input_track, input_scale)
    cascades = dgt.Digitalizer(input_cascade, input_scale)

    tracks.set_palette(0, -3, 20)
    cascades.set_palette(0, -3, 20)

    tracks.digitalize(22, 22)
    cascades.digitalize(22, 22)

    D_tracks = tracks.extracted
    D_cascades = cascades.extracted

    # data_entries = D_tracks[binnum]
    data_entries = D_cascades[binnum]

    # if binnum == 11: # only manual hard code part
    #     data_entries[14] = 0.1

    print(data_entries)

    bins_centers = np.array([0.5 * (x_bins[i] + x_bins[i+1]) for i in range(len(x_bins)-1)])
    popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[binnum, 10, 1])

    # # Generate enough x values to make the curves look smooth.
    xspace = np.linspace(2, 55, 100000)

    # # Plot the histogram and the fitted function.
    plt.bar(bins_centers, data_entries, width=x_bins[1] - x_bins[0], color='navy', label=r'Histogram entries')
    plt.plot(xspace, gaussian(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')
    plt.savefig("./MigMatPlots/Cascades/CascadeMigMatBin{}.png".format(binnum))
    plt.close()

for i in range(22):z

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

    return

# plot_zenith_errors()


# def plot_energy_reco():
#     # plt.scatter(anl.res_true, anl.res_reco)
#     heatmap, xedges, yedges = np.histogram2d(anl.res_true, anl.res_reco, bins=(x_bins, x_bins))
#     extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#     plt.xlim(2, 53)
#     plt.ylim(2, 53)
#     plt.clf()
#     plt.imshow(heatmap.T, extent=extent, origin='lower', norm = LogNorm())
#     # plt.xscale("log")
#     # plt.yscale("log")
#     plt.savefig("ICMC_with_ORCA_Reco")
#     plt.show()

def plot_energy_reco():
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(anl.res_e_true, anl.res_e_reco, bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.xscale("log")
    plt.yscale("log")
    # plt.show()
    plt.savefig("ICMC_with_ORCA_Reco")

# plot_energy_reco()

def plot_zenith_reco():
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(np.cos(anl.res_zen_true), np.cos(anl.res_zen_reco), bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()
    plt.savefig("ICMC_coszen_with_ORCA_Reco")

# plot_zenith_reco()
