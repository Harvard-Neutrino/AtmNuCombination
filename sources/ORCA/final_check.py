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

plt.style.use('./mystyle.mplstyle')

# first the plot that just reads out the digitization

def simple_readout():
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)

    IC_bins_means = np.sqrt(x[1:] * x[:-1])


    tracks = dgt.Digitalizer(input_track, input_scale)
    cascades = dgt.Digitalizer(input_cascade, input_scale)

    tracks.set_palette(0, -3, 20)
    cascades.set_palette(0, -3, 20)

    tracks.digitalize(22, 22)
    cascades.digitalize(22, 22)

    D_tracks = np.array(tracks.extracted)
    D_tracks[11][14] = 0.1 # only manual hardcode part
    D_tracks[14][3] = D_tracks[14][2] # only manual hardcode part


    D_cascades = np.array(cascades.extracted)

    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (20, 10), constrained_layout = True)
    ax1 = axes[0]
    ax2 = axes[1]
    im1 = ax1.pcolor(X, Y, D_tracks.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    im2 = ax2.pcolor(X, Y, D_cascades.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)

    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax2.set_xscale('log')
    ax2.set_yscale('log')


    fig.colorbar(im1, orientation = "vertical", format = LogFormatterMathtext(), alpha = 0.7)

    # plt.show()
    fig.savefig("./final_check_plots/simple_readout_from_digitizer")

# simple_readout()

def readout_from_fit():
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)

    bins_centers = np.sqrt(x[1:] * x[:-1])

    tracks = dgt.Digitalizer(input_track, input_scale)
    cascades = dgt.Digitalizer(input_cascade, input_scale)

    tracks.set_palette(0, -3, 20)
    cascades.set_palette(0, -3, 20)

    tracks.digitalize(22, 22)
    cascades.digitalize(22, 22)

    D_tracks = np.array(tracks.extracted)
    D_tracks[11][14] = 0.1 # only manual hardcode part
    D_tracks[14][3] = D_tracks[14][2] # only manual hardcode part


    D_cascades = np.array(cascades.extracted)

    track_data = np.zeros_like(D_tracks)
    cascade_data = np.zeros_like(D_cascades)

    track_50 = np.zeros((22,))
    track_sig = np.zeros((22,))
    track_A = np.zeros((22,))


    cascade_50 = np.zeros((22,))
    cascade_sig = np.zeros((22,))


    for i in range(22):
        t_popt, t_pcov = curve_fit(gaussian, xdata=bins_centers, ydata=D_tracks[i], p0=[i, 10, 1])
        track_50[i] = t_popt[0]
        track_sig[i] = t_popt[1]
        track_A[i] = t_popt[2]
        c_popt, c_pcov = curve_fit(gaussian, xdata=bins_centers, ydata=D_cascades[i], p0=[i, 10, 1])
        cascade_50[i] = c_popt[0]

        for j in range(22):
            track_data[i][j] = max(0.002, gaussian(bins_centers[j], *t_popt))

            cascade_data[i][j] = max(0.002, gaussian(bins_centers[j], *c_popt))

    print(track_50)
    print(track_sig)
    print(track_A)

    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (20, 10), constrained_layout = True)
    ax1 = axes[0]
    ax2 = axes[1]
    im1 = ax1.pcolor(X, Y, track_data.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    im2 = ax2.pcolor(X, Y, cascade_data.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)

    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax2.set_xscale('log')
    ax2.set_yscale('log')

    fig.colorbar(im1, orientation = "vertical", format = LogFormatterMathtext(), alpha = 0.7)

    # plt.show()
    fig.savefig("./final_check_plots/readout_from_fitter")

readout_from_fit()