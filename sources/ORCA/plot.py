import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
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

    data_entries = D_tracks[binnum]
    # data_entries = D_cascades[binnum]

    print(data_entries)

    bins_centers = np.array([0.5 * (x_bins[i] + x_bins[i+1]) for i in range(len(x_bins)-1)])
    popt, pcov = curve_fit(gaussian, xdata=bins_centers, ydata=data_entries, p0=[binnum, 20])

    # # Generate enough x values to make the curves look smooth.
    xspace = np.linspace(2, 55, 100000)

    # # Plot the histogram and the fitted function.
    plt.bar(bins_centers, data_entries, width=x_bins[1] - x_bins[0], color='navy', label=r'Histogram entries')
    plt.plot(xspace, gaussian(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')
    plt.savefig("./MigMatPlots/TrackMigMatBin{}.png".format(binnum))
    plt.close()

# for i in range(22):
#     plot_mig_hist(i)

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

# This function has some problems: heatmap definition is understood wrong I believe
# def plot_energy_reco(gaus, bins):
#     def find_reco_energy_prob(gaus, bins, true_energy, reco_energy):
#         # set bounds on IC energy
#         if true_energy > 53 or true_energy < 1.85:
#             return -1
#         # find which E_true_ORCA bin this E_true belongs to
#         for i in range(len(bins) - 1):
#             if true_energy <= bins[i + 1]:
#                 bin_num = i
#                 print("in find reco, i is ", i)
#                 break

#         # now generate a random reco number from t
#         [sigma, mu] = gaus[bin_num]
#         print("and now sigma and mu are ", sigma, mu)
#         E_reco = util.gaussian(reco_energy, mu, sigma)
#         print("and the reco is ", E_reco)

#         # return this fake ORCA MC energy
#         return E_reco
#     x = np.logspace(np.log10(2), np.log10(50), 22)
#     y = np.logspace(np.log10(2), np.log10(50), 22)
#     print("x and y are ", x, y)
#     mesh = np.zeros((22, 22), float)
#     for i in range(22):
#         for j in range(22):
#             true_energy = x[i]
#             print("currently i and j, and x and y are ", i, j, x[i], y[j])
#             res = find_reco_energy_prob(gaus, bins, x[i], y[j])
#             mesh[i][j] = res

#     fig, ax = plt.subplots()
#     ax.imshow(mesh)
#     # plt.xscale('log')
#     # plt.yscale('log')
#     plt.xlim(0.1, 50)
#     plt.ylim(0.1, 50)
#     plt.show()

# plot_energy_reco(anl.tracks.gaussians, x_bins)

def plot_energy_reco():
    # plt.scatter(anl.res_true, anl.res_reco)
    heatmap, xedges, yedges = np.histogram2d(anl.res_true, anl.res_reco, bins=(x_bins, x_bins))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.xlim(2, 53)
    plt.ylim(2, 53)
    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower', norm = LogNorm())
    # plt.xscale("log")
    # plt.yscale("log")
    plt.savefig("ICMC_with_ORCA_Reco")
    plt.show()

plot_energy_reco()
