import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit
import math
from statistics import median
from scipy import interpolate

import digitalizer as dgt
import util
from util import gaussian
from util import twogaussian
# import analyze as anl
from params import *
import NewEffective as eff

# matplotlib.rcParams.update({'font.size': 18})
plt.style.use('./paper.mplstyle')

ORCA = pd.read_csv("15x_with_interm.csv")
IC = pd.read_csv("neutrino_mc.csv")

# first plot the energy resolution
def energy_resolution():
    # first define the energy bins
    x = np.logspace(np.log10(1.85), np.log10(53), 23)
    y = np.logspace(np.log10(1.85), np.log10(53), 23)
    X, Y = np.meshgrid(x, y)
    # first Icecube
        # tracks
    IC_pid = IC["pid"]
    IC_et = IC["true_energy"]
    IC_er = IC["reco_energy"]
    IC_track_mask = IC_pid == 1
    ICtZ, xedges, yedges = np.histogram2d(IC_et[IC_track_mask], IC_er[IC_track_mask], bins=(x, y))
    for i in range(22):
        currcol = ICtZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

        # cascades
    IC_cas_mask = IC_pid == 0
    ICcZ, xedges, yedges = np.histogram2d(IC_et[IC_cas_mask], IC_er[IC_cas_mask], bins=(x, y))
    for i in range(22):
        currcol = ICcZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
        
    # ORCA
    ORCA_pid = ORCA["pid"]
    ORCA_et = ORCA["true_energy"]
    ORCA_er = ORCA["reco_energy"]
    ORCA_track_mask = ORCA_pid == 1
    ORCAtZ, xedges, yedges = np.histogram2d(ORCA_et[ORCA_track_mask], ORCA_er[ORCA_track_mask], bins=(x, y))
    for i in range(22):
        currcol = ORCAtZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

        # cascades
    ORCA_cas_mask = ORCA_pid == 0
    ORCAcZ, xedges, yedges = np.histogram2d(ORCA_et[ORCA_cas_mask], ORCA_er[ORCA_cas_mask], bins=(x, y))
    for i in range(22):
        currcol = ORCAcZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot
    
    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (20, 20), constrained_layout = True)
    ax1 = axes[0][0]
    ax1.title.set_text("IceCube Upgrade Tracks")
    ax2 = axes[0][1]
    ax2.title.set_text("IceCube Upgrade Cascades")
    im1 = ax1.pcolor(X, Y, ICtZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    im2 = ax2.pcolor(X, Y, ICcZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    ax1.set_xlim(2, 53)
    ax1.set_xlabel("True Energy [GeV]")
    ax1.set_ylim(2, 53)
    ax1.set_ylabel("Reconstructed Energy [GeV]")
    ax2.set_xlim(2, 53)
    ax2.set_xlabel("True Energy [GeV]")
    ax2.set_ylim(2, 53)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax3 = axes[1][0]
    ax3.title.set_text("ORCA Tracks")
    ax4 = axes[1][1]
    ax4.title.set_text("ORCA Cascades")
    im3 = ax3.pcolor(X, Y, ORCAtZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    im4 = ax4.pcolor(X, Y, ORCAcZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    ax3.set_xlim(2, 53)
    ax3.set_xlabel("True Energy [GeV]")
    ax3.set_ylim(2, 53)
    ax3.set_ylabel("Reconstructed Energy [GeV]")
    ax4.set_xlim(2, 53)
    ax4.set_xlabel("True Energy [GeV]")
    ax4.set_ylim(2, 53)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax4.set_xscale("log")
    ax4.set_yscale("log")

    # fig.suptitle("IceCube Upgrade and ORCA Energy Reconstruction Resolution")
    fig.colorbar(im1, ax=axes.ravel().tolist(), orientation = "vertical", format = LogFormatterMathtext(), alpha = 0.7)

    # plt.constrained_layout()
    # plt.subplots_adjust(top = 1.3)
    
    # plt.show()
    plt.savefig("./new_paper_plots/Energy_Resolution_new.png", bbox_inches = 'tight', pad_inches = 2)
    plt.close()

# energy_resolution()

def IC_track_error():
    IC_pid = IC["pid"]
    IC_et = IC["true_energy"]
    IC_er = IC["reco_energy"]
    track_mask = IC_pid == 1
    cascade_mask = IC_pid == 0

    x = np.logspace(np.log10(1.85), np.log10(53), 23)

    all_bins = []
    Cas_all_bins = []
    for i in range(22):
        all_bins.append([])
        Cas_all_bins.append([])


    median_errors = np.zeros((22,))
    top15 = np.zeros((22,))
    bottom15 = np.zeros((22,))

    Cascade_median_errors = np.zeros((22,))
    Cascade_top15 = np.zeros((22,))
    Cascade_bottom15 = np.zeros((22,))
    ICetTrack = np.array(IC_et[track_mask][:])
    ICerTrack = np.array(IC_er[track_mask][:])

    ICetCascade = np.array(IC_et[cascade_mask][:])
    ICerCascade = np.array(IC_er[cascade_mask][:])

    # print(IC_et[track_mask])
    # print(ICetTrack)


    for i in range(len(IC_er[track_mask])):
        curcol = 0
        found_place = False
        for j in range(len(x) - 1):
            # find which index we belong to
            if ICetTrack[i] < x[j + 1]:
                curcol = j 
                found_place = True
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            all_bins[curcol].append(ICerTrack[i] - ICetTrack[i])

    for i in range(len(IC_er[cascade_mask])):
        curcol = 0
        found_place = False
        for j in range(len(x) - 1):
            # find which index we belong to
            if ICetCascade[i] < x[j + 1]:
                curcol = j 
                found_place = True
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            Cas_all_bins[curcol].append(ICerCascade[i] - ICetCascade[i])

    # now select the mean from each list
    for i in range(22):
        # print(all_bins[i])
        median_errors[i] = median(all_bins[i])
        sorted_current_bin = all_bins[i].sort()
        top15[i] = abs(all_bins[i][math.ceil(0.85 * len(all_bins[i]))] - median_errors[i])
        bottom15[i] = abs(all_bins[i][math.ceil(0.15 * len(all_bins[i]))] - median_errors[i])

        Cascade_median_errors[i] = median(Cas_all_bins[i])
        Cascade_sorted_current_bin = Cas_all_bins[i].sort()
        Cascade_top15[i] = abs(Cas_all_bins[i][math.ceil(0.85 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])
        Cascade_bottom15[i] = abs(Cas_all_bins[i][math.ceil(0.15 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])

    # print(median_errors)
    # print(median_errors)
    # print(top15)
    # print(bottom15)
    IC_track_err = np.array([np.array(bottom15), np.array(top15)])
    IC_cas_err = np.array([np.array(Cascade_bottom15), np.array(Cascade_top15)])


    fig, [ax, ax2] = plt.subplots(figsize = (18, 8), ncols = 2, nrows = 1, sharey = True)
    fig.suptitle("IceCube Energy Reconstruction Errors")
    ax.set_xscale("log")
    ax.set_ylim(-25, 10)
    ax.errorbar(x[1:], median_errors, yerr = IC_track_err, fmt = 'o', capsize = 3, label = r"IceCube Track Energy Reco Median Error ($15\%$ and $85\%$)")
    ax.legend()
    ax.set_xlabel("True Energy [GeV]")
    ax2.set_xlabel("True Energy [GeV]")
    ax.set_ylabel("Energy Error [GeV]")
    # ax2.set_ylabel("Energy Error [GeV]")
    ax2.set_xscale("log")
    # ax2.set_ylim(-25, 5)
    ax2.errorbar(x[1:], Cascade_median_errors, yerr = IC_cas_err, fmt = 'o', capsize = 3, label = r"IceCube Cascade Energy Reco Median Error ($15\%$ and $85\%$)")
    plt.subplots_adjust(wspace=0, hspace=0)
    ax2.legend()
    # plt.show()
    plt.savefig("./new_paper_plots/IceCube_Energy_Errors.png", bbox_inches = 'tight', pad_inches = 1)
    plt.close()

# IC_track_error()

def ORCA_track_error():
    lo = 1.85
    hi = 53
    IC_pid = ORCA["pid"]
    IC_et = ORCA["true_energy"]
    IC_er = ORCA["reco_energy"]
    track_mask = IC_pid == 1
    cascade_mask = IC_pid == 0

    x = np.logspace(np.log10(lo), np.log10(hi), 23)

    y = np.logspace(np.log10(lo), np.log10(hi), 23)
    X, Y = np.meshgrid(x, y)
    # define the plotting start and end points
    plot_start = np.sqrt(x[0] * x[1])
    plot_end = np.sqrt(x[-1] * x[-2])

    IC_track_mask = IC_pid == 1
    ICtZ, xedges, yedges = np.histogram2d(IC_et[IC_track_mask], IC_er[IC_track_mask], bins=(x, y))
    for i in range(22):
        currcol = ICtZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

        # cascades
    IC_cas_mask = IC_pid == 0
    ICcZ, xedges, yedges = np.histogram2d(IC_et[IC_cas_mask], IC_er[IC_cas_mask], bins=(x, y))
    for i in range(22):
        currcol = ICcZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

    error_x = np.logspace(np.log10(lo), np.log10(hi), 23)

    IC_bins_means = np.sqrt(x[1:] * x[:-1])

    all_bins = []
    Cas_all_bins = []
    for i in range(22):
        all_bins.append([])
        Cas_all_bins.append([])


    median_errors = np.zeros((22,))
    top15 = np.zeros((22,))
    bottom15 = np.zeros((22,))

    Cascade_median_errors = np.zeros((22,))
    Cascade_top15 = np.zeros((22,))
    Cascade_bottom15 = np.zeros((22,))
    ICetTrack = np.array(IC_et[track_mask][:])
    ICerTrack = np.array(IC_er[track_mask][:])

    ICetCascade = np.array(IC_et[cascade_mask][:])
    ICerCascade = np.array(IC_er[cascade_mask][:])

    # print(IC_et[track_mask])
    # print(ICetTrack)


    for i in range(len(IC_er[track_mask])):
        if i <= 10:
            print(ICetTrack[i])
        curcol = 0
        found_place = False
        for j in range(len(error_x) - 1):
            # find which index we belong to
            if ICetTrack[i] < error_x[j + 1]:
                curcol = j 
                found_place = True
                if i < 10:
                    print(j)
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            all_bins[curcol].append(ICerTrack[i] - ICetTrack[i])
    # print(all_bins[0][:10])
    for i in range(len(IC_er[cascade_mask])):
        curcol = 0
        found_place = False
        for j in range(len(x) - 1):
            # find which index we belong to
            if ICetCascade[i] < x[j + 1]:
                curcol = j 
                found_place = True
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            Cas_all_bins[curcol].append(ICerCascade[i] - ICetCascade[i])

    # now select the mean from each list
    for i in range(22):
        # print(all_bins[i])
        median_errors[i] = median(all_bins[i])
        sorted_current_bin = all_bins[i].sort()
        top15[i] = abs(all_bins[i][math.ceil(0.85 * len(all_bins[i]))] - median_errors[i])
        bottom15[i] = abs(all_bins[i][math.ceil(0.15 * len(all_bins[i]))] - median_errors[i])

        Cascade_median_errors[i] = median(Cas_all_bins[i])
        Cascade_sorted_current_bin = Cas_all_bins[i].sort()
        Cascade_top15[i] = abs(Cas_all_bins[i][math.ceil(0.85 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])
        Cascade_bottom15[i] = abs(Cas_all_bins[i][math.ceil(0.15 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])

    # print(median_errors)
    # print(median_errors)
    # print(top15)
    # print(bottom15)
    IC_track_err = np.array([np.array(bottom15), np.array(top15)])
    IC_cas_err = np.array([np.array(Cascade_bottom15), np.array(Cascade_top15)])


    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 10), constrained_layout = True, sharex = True, sharey = 'row', gridspec_kw = {'height_ratios':[1, 2]})
    # fig, axes = plt.subplots(figsize = (15, 10), ncols = 2, nrows = 2, sharey = True)
    fig.suptitle("ORCA Energy Reconstruction Errors")
    ax, ax2 = axes[0][0], axes[0][1]
    ax3, ax4 = axes[1][0], axes[1][1]
    ax.set_xscale("log")
    # ax.set_ylim(-25, 10)
    ax.errorbar(IC_bins_means, median_errors, yerr = IC_track_err, fmt = 'o', capsize = 3, label = r"ORCA Track Energy Reco Median Error ($15\%$ and $85\%$)")
    ax.legend(fontsize = 13)
    ax.set_xlabel("True Energy [GeV]")
    ax2.set_xlabel("True Energy [GeV]")
    ax.set_ylabel("Energy Error [GeV]")
    # ax2.set_ylabel("Energy Error [GeV]")
    ax2.set_xscale("log")
    # ax2.set_ylim(-25, 5)
    ax2.errorbar(IC_bins_means, Cascade_median_errors, yerr = IC_cas_err, fmt = 'o', capsize = 3, label = r"ORCA Cascade Energy Reco Median Error ($15\%$ and $85\%$)")
    # plt.subplots_adjust(wspace=0, hspace=0)
    ax2.legend(fontsize = 13)

    im1 = ax3.pcolor(X, Y, ICtZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    im2 = ax4.pcolor(X, Y, ICcZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    ax3.set_xlim(plot_start, plot_end)
    ax3.set_xlabel("True Energy [GeV]")
    ax3.set_ylim(lo, hi)
    ax3.set_ylabel("Reconstructed Energy [GeV]")
    ax4.set_xlim(plot_start, plot_end)
    ax4.set_xlabel("True Energy [GeV]")
    ax4.set_ylim(lo, hi)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    fig.colorbar(im1, orientation = "vertical", format = LogFormatterMathtext(), alpha = 0.7)

    # plt.show()

    plt.savefig("./new_paper_plots/ORCA_Energy_Errors.png", bbox_inches = 'tight', pad_inches = 1)
    plt.close()

# ORCA_track_error()

def IC_resolution_with_errors():
    lo = 0.5
    hi = 105
    # first define the energy bins
    x = np.logspace(np.log10(lo), np.log10(hi), 23)
    y = np.logspace(np.log10(lo), np.log10(hi), 23)
    X, Y = np.meshgrid(x, y)
    # define the plotting start and end points
    plot_start = np.sqrt(x[0] * x[1])
    plot_end = np.sqrt(x[-1] * x[-2])
    # print(plot_start, plot_end)

    # Here are the IC resolution plots
        # tracks
    IC_pid = IC["pid"]
    IC_et = IC["true_energy"]
    IC_er = IC["reco_energy"]
    IC_track_mask = IC_pid == 1
    ICtZ, xedges, yedges = np.histogram2d(IC_et[IC_track_mask], IC_er[IC_track_mask], bins=(x, y))
    for i in range(22):
        currcol = ICtZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

        # cascades
    IC_cas_mask = IC_pid == 0
    ICcZ, xedges, yedges = np.histogram2d(IC_et[IC_cas_mask], IC_er[IC_cas_mask], bins=(x, y))
    for i in range(22):
        currcol = ICcZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

    # Here are the IC error bars
    IC_pid = IC["pid"]
    IC_et = IC["true_energy"]
    IC_er = IC["reco_energy"]
    track_mask = IC_pid == 1
    cascade_mask = IC_pid == 0

    # THIS ENDS THE RESOLUTION PLOT

    error_x = np.logspace(np.log10(lo), np.log10(hi), 23)
    IC_bins_means = np.sqrt(x[1:] * x[:-1])

    all_bins = []
    Cas_all_bins = []
    for i in range(22):
        all_bins.append([])
        Cas_all_bins.append([])


    median_errors = np.zeros((22,))
    top15 = np.zeros((22,))
    bottom15 = np.zeros((22,))
    track_median_ratio = np.zeros((22,))
    track_ratio_top = np.zeros((22,))
    track_ratio_bottom = np.zeros((22,))

    Cascade_median_errors = np.zeros((22,))
    cascade_median_ratio = np.zeros((22,))
    Cascade_top15 = np.zeros((22,))
    Cascade_bottom15 = np.zeros((22,))
    cascade_ratio_top = np.zeros((22,))
    cascade_ratio_bottom = np.zeros((22,))


    ICetTrack = np.array(IC_et[track_mask][:])
    ICerTrack = np.array(IC_er[track_mask][:])

    ICetCascade = np.array(IC_et[cascade_mask][:])
    ICerCascade = np.array(IC_er[cascade_mask][:])

    for i in range(len(IC_er[track_mask])):
        found_place = False
        curcol = 0
        for j in range(len(error_x) - 1):
            # find which index we belong to
            if ICetTrack[i] < error_x[j + 1]:
                curcol = j 
                found_place = True
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            all_bins[curcol].append(ICerTrack[i] - ICetTrack[i])

    for i in range(len(IC_er[cascade_mask])):
        curcol = 0
        found_place = False
        for j in range(len(x) - 1):
            # find which index we belong to
            if ICetCascade[i] < error_x[j + 1]:
                curcol = j 
                found_place = True
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            Cas_all_bins[curcol].append(ICerCascade[i] - ICetCascade[i])

    # print(len(all_bins[0]))
    # print(len(Cas_all_bins[0]))

    # print(len(all_bins[1]))
    # print(len(Cas_all_bins[1]))

    # print(all_bins[0])
    # print(Cas_all_bins[0])

    # now select the mean from each list
    for i in range(22):
        # print(all_bins[i])
        median_errors[i] = median(all_bins[i])
        sorted_current_bin = all_bins[i].sort()
        top15[i] = abs(all_bins[i][math.floor(0.85 * len(all_bins[i]))] - median_errors[i])
        bottom15[i] = abs(all_bins[i][math.ceil(0.15 * len(all_bins[i]))] - median_errors[i])

        track_ratio_top[i] = top15[i] / IC_bins_means[i]
        track_ratio_bottom[i] = bottom15[i] / IC_bins_means[i]
        track_median_ratio[i] = median_errors[i] / IC_bins_means[i]

        Cascade_median_errors[i] = median(Cas_all_bins[i])
        Cascade_sorted_current_bin = Cas_all_bins[i].sort()
        Cascade_top15[i] = abs(Cas_all_bins[i][math.floor(0.85 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])
        Cascade_bottom15[i] = abs(Cas_all_bins[i][math.ceil(0.15 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])

        cascade_ratio_top[i] = Cascade_top15[i] / IC_bins_means[i]
        cascade_ratio_bottom[i] = Cascade_bottom15[i] / IC_bins_means[i]
        cascade_median_ratio[i] = Cascade_median_errors[i] / IC_bins_means[i]


    IC_track_err = np.array([np.array(bottom15), np.array(top15)])
    IC_cas_err = np.array([np.array(Cascade_bottom15), np.array(Cascade_top15)])

    IC_track_ratio_err = np.array([np.array(track_ratio_bottom), np.array(track_ratio_top)])
    IC_cas_ratio_err = np.array([np.array(cascade_ratio_bottom), np.array(cascade_ratio_top)])

    # here process the lists of line points that go into the resolutioon plot
    
    track_median = np.zeros_like(IC_bins_means)
    track_up = np.zeros_like(IC_bins_means)
    track_down = np.zeros_like(IC_bins_means)
    for i in range(len(IC_bins_means)):
        track_median[i] = IC_bins_means[i] + median_errors[i]
        track_up[i] = IC_bins_means[i] + median_errors[i]+ top15[i]
        track_down[i] = IC_bins_means[i] + median_errors[i]- bottom15[i]

    cascade_median = np.zeros_like(IC_bins_means)
    cascade_up = np.zeros_like(IC_bins_means)
    cascade_down = np.zeros_like(IC_bins_means)
    for i in range(len(IC_bins_means)):
        cascade_median[i] = IC_bins_means[i] + Cascade_median_errors[i]
        cascade_up[i] = IC_bins_means[i] + Cascade_median_errors[i]+ Cascade_top15[i]
        cascade_down[i] = IC_bins_means[i] + Cascade_median_errors[i]- Cascade_bottom15[i]

    # Now interpolate this line in the histogram
    # print(IC_bins_means)
    f_track_median = interpolate.interp1d(IC_bins_means, track_median, kind = "quadratic")
    f_track_top = interpolate.interp1d(IC_bins_means, track_up, kind = "quadratic")
    f_track_bot = interpolate.interp1d(IC_bins_means, track_down, kind = "quadratic")
    f_cascade_median = interpolate.interp1d(IC_bins_means, cascade_median, kind = "quadratic")
    f_cascade_top = interpolate.interp1d(IC_bins_means, cascade_up, kind = "quadratic")
    f_cascade_bot = interpolate.interp1d(IC_bins_means, cascade_down, kind = "quadratic")

    plot_error_x = np.logspace(np.log10(IC_bins_means[0]), np.log10(IC_bins_means[-1]), 1000)
    plot_track_median = f_track_median(plot_error_x)
    plot_track_top = f_track_top(plot_error_x)
    plot_track_bot = f_track_bot(plot_error_x)
    plot_cascade_median = f_cascade_median(plot_error_x)
    plot_cascade_top = f_cascade_top(plot_error_x)
    plot_cascade_bot = f_cascade_bot(plot_error_x)

    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (12, 8), constrained_layout = True, sharex = True, sharey = 'row', gridspec_kw = {'height_ratios':[1, 2]})
    ax3 = axes[0][0]
    ax4 = axes[0][1]
    ax1 = axes[1][0]
    ax2 = axes[1][1]
    im1 = ax1.pcolor(X, Y, ICtZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    im2 = ax2.pcolor(X, Y, ICcZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    ax1.set_xlim(plot_start, plot_end)
    ax1.set_xlabel("True Energy [GeV]")
    ax1.set_ylim(lo, hi)
    ax1.set_ylabel("Reconstructed Energy [GeV]")
    ax2.set_xlim(plot_start, plot_end)
    ax2.set_xlabel("True Energy [GeV]")
    ax2.set_ylim(lo, hi)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    fig.colorbar(im1, orientation = "vertical", format = LogFormatterMathtext(), alpha = 0.7)

    # also plot the error bars
    ax1.plot(plot_error_x, plot_track_median, label = r"\textbf{Median} 50\%", color = 'navy', alpha = 0.7)
    ax1.plot(plot_error_x, plot_track_top, label = r"\textbf{Top} 15\%", color = 'navy', linestyle = '--', alpha = 0.7)
    ax1.plot(plot_error_x, plot_track_bot, label = r"\textbf{Bottom} 15\%", color = 'navy', linestyle = '-.', alpha = 0.7)
    ax1.legend(fontsize = 14)


    ax2.plot(plot_error_x, plot_cascade_median, label = r"\textbf{Median} 50\%", color = 'navy', alpha = 0.7)
    ax2.plot(plot_error_x, plot_cascade_top, label = r"\textbf{Top} 15\%", color = 'navy', linestyle = '--', alpha = 0.7)
    ax2.plot(plot_error_x, plot_cascade_bot, label = r"\textbf{Bottom} 15\%", color = 'navy', linestyle = '-.', alpha = 0.7)
    ax2.legend(fontsize = 14)

    # ax3.set_ylim(-25, 10)
    # ax3.errorbar(IC_bins_means, track_median_ratio, yerr = IC_track_ratio_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax3.errorbar(IC_bins_means, median_errors, yerr = IC_track_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax3.legend(fontsize = 14)
    # ax3.set_xlabel("True Energy [GeV]")
    # ax4.set_xlabel("True Energy [GeV]")
    ax3.set_ylabel("Energy Error[GeV]")
    # ax4.set_xscale("log")
    # ax2.set_ylim(-25, 5)
    # ax4.errorbar(IC_bins_means, cascade_median_ratio, yerr = IC_cas_ratio_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax4.errorbar(IC_bins_means, Cascade_median_errors, yerr = IC_cas_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    # plt.subplots_adjust(hspace=0)
    ax4.legend(fontsize = 14)
    # plt.subplots_adjust(wspace=0.1, hspace=0.03)

    # plt.show()
    plt.savefig("./new_paper_plots/IC_energy_Reconstruction", bbox_inches = 'tight', pad_inches = 0.3)

# IC_resolution_with_errors()

def ORCA_resolution_with_errors():
    lo = 1.85
    hi = 53
    # first define the energy bins
    x = np.logspace(np.log10(lo), np.log10(hi), 23)
    y = np.logspace(np.log10(lo), np.log10(hi), 23)
    X, Y = np.meshgrid(x, y)
    # define the plotting start and end points
    plot_start = np.sqrt(x[0] * x[1])
    plot_end = np.sqrt(x[-1] * x[-2])
    # print(plot_start, plot_end)

    # Here are the IC resolution plots
        # tracks
    IC_pid = ORCA["pid"]
    IC_et = ORCA["true_energy"]
    IC_er = ORCA["reco_energy"]
    IC_track_mask = IC_pid == 1
    ICtZ, xedges, yedges = np.histogram2d(IC_et[IC_track_mask], IC_er[IC_track_mask], bins=(x, y))
    for i in range(22):
        currcol = ICtZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

        # cascades
    IC_cas_mask = IC_pid == 0
    ICcZ, xedges, yedges = np.histogram2d(IC_et[IC_cas_mask], IC_er[IC_cas_mask], bins=(x, y))
    for i in range(22):
        currcol = ICcZ[i][:]
        tot = 0
        for j in range(22):
            tot += currcol[j]
        for j in range(22):
            currcol[j] = currcol[j] / tot

    # Here are the IC error bars
    IC_pid = ORCA["pid"]
    IC_et = ORCA["true_energy"]
    IC_er = ORCA["reco_energy"]
    track_mask = IC_pid == 1
    cascade_mask = IC_pid == 0

    # THIS ENDS THE RESOLUTION PLOT

    error_x = np.logspace(np.log10(lo), np.log10(hi), 23)
    IC_bins_means = np.sqrt(x[1:] * x[:-1])

    all_bins = []
    Cas_all_bins = []
    for i in range(22):
        all_bins.append([])
        Cas_all_bins.append([])


    median_errors = np.zeros((22,))
    top15 = np.zeros((22,))
    bottom15 = np.zeros((22,))
    # track_median_ratio = np.zeros((22,))
    # track_ratio_top = np.zeros((22,))
    # track_ratio_bottom = np.zeros((22,))

    Cascade_median_errors = np.zeros((22,))
    # cascade_median_ratio = np.zeros((22,))
    Cascade_top15 = np.zeros((22,))
    Cascade_bottom15 = np.zeros((22,))
    # cascade_ratio_top = np.zeros((22,))
    # cascade_ratio_bottom = np.zeros((22,))


    ICetTrack = np.array(IC_et[track_mask][:])
    ICerTrack = np.array(IC_er[track_mask][:])

    ICetCascade = np.array(IC_et[cascade_mask][:])
    ICerCascade = np.array(IC_er[cascade_mask][:])

    for i in range(len(IC_er[track_mask])):
        # if i <= 10:
        #     print(ICetTrack[i])
        found_place = False
        curcol = 0
        for j in range(len(error_x) - 1):
            # find which index we belong to
            if ICetTrack[i] < error_x[j + 1]:
                curcol = j 
                found_place = True
                # if i < 10:
                #     print(j)
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            all_bins[curcol].append(ICerTrack[i] - ICetTrack[i])
    # print(all_bins[0][:10])
    for i in range(len(IC_er[cascade_mask])):
        curcol = 0
        found_place = False
        for j in range(len(x) - 1):
            # find which index we belong to
            if ICetCascade[i] < error_x[j + 1]:
                curcol = j 
                found_place = True
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            Cas_all_bins[curcol].append(ICerCascade[i] - ICetCascade[i])

    # print(len(all_bins[0]))
    # print(len(Cas_all_bins[0]))

    # print(len(all_bins[1]))
    # print(len(Cas_all_bins[1]))

    # print(all_bins[0])
    # print(Cas_all_bins[0])

    # now select the mean from each list
    for i in range(22):
        # print(all_bins[i])
        median_errors[i] = median(all_bins[i])
        sorted_current_bin = all_bins[i].sort()
        top15[i] = abs(all_bins[i][math.floor(0.85 * len(all_bins[i]))] - median_errors[i])
        bottom15[i] = abs(all_bins[i][math.ceil(0.15 * len(all_bins[i]))] - median_errors[i])

        # track_ratio_top[i] = top15[i] / IC_bins_means[i]
        # track_ratio_bottom[i] = bottom15[i] / IC_bins_means[i]
        # track_median_ratio[i] = median_errors[i] / IC_bins_means[i]

        Cascade_median_errors[i] = median(Cas_all_bins[i])
        Cascade_sorted_current_bin = Cas_all_bins[i].sort()
        Cascade_top15[i] = abs(Cas_all_bins[i][math.floor(0.85 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])
        Cascade_bottom15[i] = abs(Cas_all_bins[i][math.ceil(0.15 * len(Cas_all_bins[i]))] - Cascade_median_errors[i])

        # cascade_ratio_top[i] = Cascade_top15[i] / IC_bins_means[i]
        # cascade_ratio_bottom[i] = Cascade_bottom15[i] / IC_bins_means[i]
        # cascade_median_ratio[i] = Cascade_median_errors[i] / IC_bins_means[i]

    # print(median_errors)
    # print(Cascade_median_errors)

    IC_track_err = np.array([np.array(bottom15), np.array(top15)])
    IC_cas_err = np.array([np.array(Cascade_bottom15), np.array(Cascade_top15)])

    # IC_track_ratio_err = np.array([np.array(track_ratio_bottom), np.array(track_ratio_top)])
    # IC_cas_ratio_err = np.array([np.array(cascade_ratio_bottom), np.array(cascade_ratio_top)])

    # here process the lists of line points that go into the resolutioon plot
    
    track_median = np.zeros_like(IC_bins_means)
    track_up = np.zeros_like(IC_bins_means)
    track_down = np.zeros_like(IC_bins_means)
    for i in range(len(IC_bins_means)):
        track_median[i] = IC_bins_means[i] + median_errors[i]
        track_up[i] = IC_bins_means[i] + median_errors[i] + top15[i]
        track_down[i] = IC_bins_means[i] + median_errors[i] - bottom15[i]

    # print(track_median)
    # print(track_up)
    # print(track_down)

    cascade_median = np.zeros_like(IC_bins_means)
    cascade_up = np.zeros_like(IC_bins_means)
    cascade_down = np.zeros_like(IC_bins_means)
    for i in range(len(IC_bins_means)):
        cascade_median[i] = IC_bins_means[i] + Cascade_median_errors[i]
        cascade_up[i] = IC_bins_means[i] + Cascade_median_errors[i] + Cascade_top15[i]
        cascade_down[i] = IC_bins_means[i] + Cascade_median_errors[i] - Cascade_bottom15[i]

    # Now interpolate this line in the histogram
    # print(IC_bins_means)
    f_track_median = interpolate.interp1d(IC_bins_means, track_median, kind = "quadratic")
    f_track_top = interpolate.interp1d(IC_bins_means, track_up, kind = "quadratic")
    f_track_bot = interpolate.interp1d(IC_bins_means, track_down, kind = "quadratic")
    f_cascade_median = interpolate.interp1d(IC_bins_means, cascade_median, kind = "quadratic")
    f_cascade_top = interpolate.interp1d(IC_bins_means, cascade_up, kind = "quadratic")
    f_cascade_bot = interpolate.interp1d(IC_bins_means, cascade_down, kind = "quadratic")

    plot_error_x = np.logspace(np.log10(IC_bins_means[0]), np.log10(IC_bins_means[-1]), 1000)
    plot_track_median = f_track_median(plot_error_x)
    plot_track_top = f_track_top(plot_error_x)
    plot_track_bot = f_track_bot(plot_error_x)
    plot_cascade_median = f_cascade_median(plot_error_x)
    plot_cascade_top = f_cascade_top(plot_error_x)
    plot_cascade_bot = f_cascade_bot(plot_error_x)

    fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (12, 8), constrained_layout = True, sharex = True, sharey = 'row', gridspec_kw = {'height_ratios':[1, 2]})
    ax3 = axes[0][0]
    ax4 = axes[0][1]
    ax1 = axes[1][0]
    ax2 = axes[1][1]
    im1 = ax1.pcolor(X, Y, ICtZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    im2 = ax2.pcolor(X, Y, ICcZ.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    ax1.set_xlim(plot_start, plot_end)
    ax1.set_xlabel("True Energy [GeV]")
    ax1.set_ylim(lo, hi)
    ax1.set_ylabel("Reconstructed Energy [GeV]")
    ax2.set_xlim(plot_start, plot_end)
    ax2.set_xlabel("True Energy [GeV]")
    ax2.set_ylim(lo, hi)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    fig.colorbar(im1, orientation = "vertical", format = LogFormatterMathtext(), alpha = 0.7)

    # also plot the error bars
    ax1.plot(plot_error_x, plot_track_median, label = r"\textbf{Median} 50\%", color = 'navy', alpha = 0.7)
    ax1.plot(plot_error_x, plot_track_top, label = r"\textbf{Top} 15\%", color = 'navy', linestyle = '--', alpha = 0.7)
    ax1.plot(plot_error_x, plot_track_bot, label = r"\textbf{Bottom} 15\%", color = 'navy', linestyle = '-.', alpha = 0.7)
    ax1.legend(fontsize = 14)


    ax2.plot(plot_error_x, plot_cascade_median, label = r"\textbf{Median} 50\%", color = 'navy', alpha = 0.7)
    ax2.plot(plot_error_x, plot_cascade_top, label = r"\textbf{Top} 15\%", color = 'navy', linestyle = '--', alpha = 0.7)
    ax2.plot(plot_error_x, plot_cascade_bot, label = r"\textbf{Bottom} 15\%", color = 'navy', linestyle = '-.', alpha = 0.7)
    ax2.legend(fontsize = 14)

    # ax3.set_ylim(-25, 10)
    # ax3.errorbar(IC_bins_means, track_median_ratio, yerr = IC_track_ratio_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax3.errorbar(IC_bins_means, median_errors, yerr = IC_track_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax3.legend(fontsize = 14)
    # ax3.set_xlabel("True Energy [GeV]")
    # ax4.set_xlabel("True Energy [GeV]")
    ax3.set_ylabel("Energy Error[GeV]")
    # ax4.set_xscale("log")
    # ax2.set_ylim(-25, 5)
    # ax4.errorbar(IC_bins_means, cascade_median_ratio, yerr = IC_cas_ratio_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax4.errorbar(IC_bins_means, Cascade_median_errors, yerr = IC_cas_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    # plt.subplots_adjust(hspace=0)
    ax4.legend(fontsize = 14)
    # plt.subplots_adjust(wspace=0.08, hspace=0.03)

    # plt.show()
    plt.savefig("./new_paper_plots/ORCA_energy_Reconstruction", bbox_inches = 'tight', pad_inches = 0.3)

# ORCA_resolution_with_errors()

def ORCA_zenith_resolution():
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(np.cos(ORCA["true_zenith"]), np.cos(ORCA["reco_zenith"]), bins=(x, y))
    
    zen_t = ORCA["true_zenith"]
    zen_r = ORCA["reco_zenith"]

    plot_start = .5 * (x[0] + x[1])
    plot_end = .5 * (x[-1] + x[-2])

    error_zen = np.linspace(-1, 1, 20)
    ORCA_bins_means = 0.5 * (x[1:] + x[:-1])

    all_bins = []
    for i in range(19):
        all_bins.append([])

    median_errors = np.zeros((19,))
    top15 = np.zeros((19,))
    bottom15 = np.zeros((19,))

    zen_t = np.array(zen_t[:])
    zen_r = np.array(zen_r[:])

    for i in range(len(zen_t)):
        found_place = False
        curcol = 0
        for j in range(len(error_zen) - 1):
            # find which index we belong to
            if np.cos(zen_t[i]) < error_zen[j + 1]:
                curcol = j 
                found_place = True
                # if i < 10:
                #     print(j)
                break
            else:
                curcol = j
                continue 
        # now put the difference between reco and true into the corresponding list
        if found_place:
            all_bins[curcol].append(np.cos(zen_r[i]) - np.cos(zen_t[i]))
        
    for i in range(19):
    # print(all_bins[i])
        median_errors[i] = median(all_bins[i])
        sorted_current_bin = all_bins[i].sort()
        top15[i] = abs(all_bins[i][math.floor(0.85 * len(all_bins[i]))] - median_errors[i])
        bottom15[i] = abs(all_bins[i][math.ceil(0.15 * len(all_bins[i]))] - median_errors[i])
    
    zen_error = np.array([np.array(bottom15), np.array(top15)])

    zen_median = np.zeros_like(ORCA_bins_means)
    zen_up = np.zeros_like(ORCA_bins_means)
    zen_down = np.zeros_like(ORCA_bins_means)
    for i in range(len(ORCA_bins_means)):
        zen_median[i] = ORCA_bins_means[i] + median_errors[i]
        zen_up[i] = ORCA_bins_means[i] + median_errors[i] + top15[i]
        zen_down[i] = ORCA_bins_means[i] + median_errors[i] - bottom15[i]

    
    f_zen_median = interpolate.interp1d(ORCA_bins_means, zen_median, kind = "quadratic")
    f_zen_top = interpolate.interp1d(ORCA_bins_means, zen_up, kind = "quadratic")
    f_zen_bot = interpolate.interp1d(ORCA_bins_means, zen_down, kind = "quadratic")

    plot_error_zen = np.linspace(ORCA_bins_means[0], ORCA_bins_means[-1], 1000)
    plot_zen_median = f_zen_median(plot_error_zen)
    plot_zen_top = f_zen_top(plot_error_zen)
    plot_zen_bot = f_zen_bot(plot_error_zen)

    fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize = (6, 7.5), constrained_layout = True, sharex = True, gridspec_kw = {'height_ratios':[1, 2]})
    ax3 = axes[0]
    ax3.set_ylim(-0.2, 0.2)
    ax1 = axes[1]
    im1 = ax1.pcolor(X, Y, Z.T, cmap = "plasma", norm = LogNorm(), alpha = 0.7)
    ax1.set_xlim(plot_start, plot_end)
    ax1.set_xlabel(r"True $\cos \theta$")
    ax1.set_ylim(-1, 1)
    ax1.set_ylabel(r"Reconstructed $\cos \theta$")
    fig.colorbar(im1, orientation = "vertical", format = LogFormatterMathtext(), alpha = 0.7)

    # also plot the error bars
    ax1.plot(plot_error_zen, plot_zen_median, label = r"\textbf{Median} $50\%$", color = 'navy', alpha = 0.7)
    ax1.plot(plot_error_zen, plot_zen_top, label = r"\textbf{Top} $15\%$", color = 'navy', linestyle = '--', alpha = 0.7)
    ax1.plot(plot_error_zen, plot_zen_bot, label = r"\textbf{Bottom} $15\%$", color = 'navy', linestyle = '-.', alpha = 0.7)
    ax1.legend(fontsize = 14)

    # ax3.set_ylim(-25, 10)
    # ax3.errorbar(IC_bins_means, track_median_ratio, yerr = IC_track_ratio_err, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax3.errorbar(ORCA_bins_means, median_errors, yerr = zen_error, fmt = 'o', fillstyle = 'none', capsize = 5, label = r"Energy Reco Error (Median, $15\%$ and $85\%$)")
    ax3.legend(fontsize = 14)
    ax3.set_ylabel(r"$\cos \theta$ Error")

    # plt.show()
    plt.savefig("./new_paper_plots/Zenith_Resolution.png", bbox_inches = 'tight', pad_inches = 0.3)

# ORCA_zenith_resolution()

def range_zenith_resolution():
    zen_true = ORCA["true_zenith"]
    zen_reco = ORCA["reco_zenith"]
    e_true = ORCA["true_energy"]
    e_reco = ORCA["reco_energy"]
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    def rangezen(elo, ehi):
        select_true = []
        select_reco = []
        for i in range(len(zen_true)):
            if e_true[i] >= elo and e_true[i] <= ehi:
                select_true.append(zen_true[i])
                select_reco.append(zen_reco[i])
        Z, xedges, yedges = np.histogram2d(np.cos(select_true), np.cos(select_reco), bins=(x, y))
        return Z
    Z1 = rangezen(1, 5)
    Z2 = rangezen(5, 10)
    Z3 = rangezen(10, 50)
    fig, axes = plt.subplots(nrows = 1, ncols = 3, figsize = (35, 13), constrained_layout = True)
    ax1, ax2, ax3 = axes[0], axes[1], axes[2]
    im1 = ax1.pcolor(X, Y, Z1.T, cmap = "plasma", norm = LogNorm())
    im2 = ax2.pcolor(X, Y, Z2.T, cmap = "plasma", norm = LogNorm())
    im3 = ax3.pcolor(X, Y, Z3.T, cmap = "plasma", norm = LogNorm())
    ax1.set_xlim(-1, 1)
    ax2.set_xlim(-1, 1)
    ax3.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)
    ax2.set_ylim(-1, 1)
    ax3.set_ylim(-1, 1)
    ax1.set_title(r"$E_{{true}} \in [1, 5]$ GeV")
    ax2.set_title(r"$E_{{true}} \in [5, 10]$ GeV")
    ax3.set_title(r"$E_{{true}} \in [10, 50]$ GeV")
    ax2.set_ylabel("Reconstructed Cosine Zenith")
    ax3.set_xlabel("True Cosine Zenith")
    plt.colorbar(im1, orientation = "vertical", format = LogFormatterMathtext())
    # plt.show()
    plt.savefig("./new_paper_plots/Ranged_Zenith_Resolution.png", bbox_inches = 'tight', pad_inches = 2)

# range_zenith_resolution()


# then plot Effective Area comparisons
def effective_volumes():
    IC = eff.ICEffectiveAnalysis(1, 50, 21)
    ORCA = eff.ORCAEffectiveAnalysis(num_bins = 20)
    IC.compute()
    ORCA.compute()
    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (18, 6))
    ax1, ax2 = axes[0], axes[1]

    # first plot IC
    ehist, ebin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[0][0])
    ebarhist, ebarbin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[1][0])
    muhist, mubin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[0][1])
    mubarhist, mubarbin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[1][1])
    tauhist, taubin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[0][2])
    taubarhist, taubarbin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[1][2])
    nchist, ncbin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[0][3])
    ncbarhist, ncbarbin = np.histogram(IC.energy, bins = IC.bins, weights = IC.volumes[1][3])

    ax1.hist(ebin[:-1], ebin, weights = ehist / IC.widths, color = 'blue', label=r"$\nu_eCC$", histtype="step")
    ax1.hist(ebarbin[:-1], ebarbin, weights = ebarhist / IC.widths, color = 'blue', label=r"$\overline{\nu}_eCC$", \
                            linestyle = '--', histtype="step")
    ax1.hist(mubin[:-1], mubin, weights = muhist / IC.widths, color = 'red', label=r"$\nu_{\mu}CC$", histtype="step")
    ax1.hist(mubarbin[:-1], mubarbin, weights = mubarhist / IC.widths,  color = 'red',label=r"$\overline{\nu}_{\mu}CC$", \
                            linestyle = '--', histtype="step")
    ax1.hist(taubin[:-1], taubin, weights = tauhist / IC.widths, color = 'green', label=r"$\nu_{\tau}CC$", histtype="step")
    ax1.hist(taubarbin[:-1], taubarbin, weights = taubarhist / IC.widths, color = 'green', label=r"$\overline{\nu}_{\tau}CC$",\
                            linestyle = "--", histtype="step")
    ax1.hist(ncbin[:-1], ncbin, weights = nchist / IC.widths, color = "brown", label=r"$\nu NC$", histtype="step")
    ax1.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist / IC.widths,color = "brown", label=r"$\overline{\nu} NC$",\
                            linestyle = "--", histtype="step")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
    ax1.set_xlim(IC.lo, IC.hi)
    ax1.set_ylabel(r"Effective Volume [m$^3$]")
    ax1.grid(True)
    ax1.legend(loc = 2, fontsize = 16)
    ax1.title.set_text("IceCube Upgrade Deepcore")

    # then plot ORCA
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.evol, color = 'blue', label=r"$\nu_eCC$", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.ebarvol , color = 'blue', label=r"$\overline{\nu}_eCC$", \
                            linestyle = '--', histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.muvol , color = 'red',label=r"$\nu_{\mu}CC$", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.mubarvol , color = 'red',label=r"$\overline{\nu}_{\mu}CC$", \
                            linestyle = '--', histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.tauvol , color = 'green',label=r"$\nu_{\tau}CC$", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.taubarvol , color = 'green',label=r"$\overline{\nu}_{\tau}CC$",\
                             linestyle = "--", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.ncvol ,color = "brown", label=r"$\nu NC$", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.ncbarvol ,color = "brown", label=r"$\overline{\nu} NC$", \
                            linestyle = "--", histtype="step")
    ax2.set_xscale("log")
    # ax.set_yscale("log")
    ax2.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
    ax2.set_ylabel(r"Effective Volume [m$^3$]")
    ax2.set_xlim(1, 50)
    ax2.grid(True)
    ax2.legend(loc = 2, fontsize = 16)
    # ax1.tick_params(axis='both', which='major', labelsize=18)
    # ax2.tick_params(axis='both', which='major', labelsize=18)

    plt.subplots_adjust(wspace=0.3)
    ax2.title.set_text("ORCA")

    # fig.suptitle("IceCube Upgrade DeepCore and ORCA Effective Volumes")
    # plt.show()
    plt.savefig("./new_paper_plots/Effective_Volumes.png", bbox_inches = 'tight', pad_inches = 0.2)
    plt.close()

effective_volumes() 