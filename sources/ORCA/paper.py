import matplotlib.pyplot as plt
import matplotlib
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
import NewEffective as eff

matplotlib.rcParams.update({'font.size': 20})

ORCA = pd.read_csv(savename)
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
    im1 = ax1.pcolor(X, Y, ICtZ.T, cmap = "gray_r", norm = LogNorm())
    im2 = ax2.pcolor(X, Y, ICcZ.T, cmap = "gray_r", norm = LogNorm())
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
    im3 = ax3.pcolor(X, Y, ORCAtZ.T, cmap = "gray_r", norm = LogNorm())
    im4 = ax4.pcolor(X, Y, ORCAcZ.T, cmap = "gray_r", norm = LogNorm())
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

    fig.suptitle("IceCube Upgrade and ORCA Energy Reconstruction Resolution")
    fig.colorbar(im1, ax=axes.ravel().tolist(), orientation = "vertical", format = LogFormatterMathtext())

    # plt.constrained_layout()
    # plt.subplots_adjust(top = 1.3)
    
    # plt.show()
    plt.savefig("./paper_plots/Energy_Resolution.png")
    plt.close()

# energy_resolution()

def zenith_resolution():
    x = np.linspace(-1, 1, 20)
    y = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(np.cos(ORCA["true_zenith"]), np.cos(ORCA["reco_zenith"]), bins=(x, y))
    im = plt.pcolor(X, Y, Z.T, cmap = "gray_r", norm = LogNorm())
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.xlabel("Cosine True Zenith")
    plt.ylabel("Cosine Reco Zenith")
    plt.colorbar(im, orientation = "vertical", format = LogFormatterMathtext())
    plt.savefig("./paper_plots/Zenith_Resolution.png")
    plt.close()

zenith_resolution()

# then plot Effective Area comparisons
def effective_volumes():
    IC = eff.ICEffectiveAnalysis(1, 50, 21)
    ORCA = eff.ORCAEffectiveAnalysis(num_bins = 20)
    IC.compute()
    ORCA.compute()
    fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (25, 10))
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

    ax1.hist(ebin[:-1], ebin, weights = ehist / IC.widths, label=r"$\nu_eCC$", color="red", histtype="step")
    ax1.hist(ebarbin[:-1], ebarbin, weights = ebarhist / IC.widths, label=r"$\overline{\nu}_eCC$", \
                            linestyle = '--', color="red", histtype="step")
    ax1.hist(mubin[:-1], mubin, weights = muhist / IC.widths, label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
    ax1.hist(mubarbin[:-1], mubarbin, weights = mubarhist / IC.widths, label=r"$\overline{\nu}_{\mu}CC$", \
                            linestyle = '--', color="blue", histtype="step")
    ax1.hist(taubin[:-1], taubin, weights = tauhist / IC.widths, label=r"$\nu_{\tau}CC$", color="green", histtype="step")
    ax1.hist(taubarbin[:-1], taubarbin, weights = taubarhist / IC.widths, label=r"$\overline{\nu}_{\tau}CC$",\
                            color="green", linestyle = "--", histtype="step")
    ax1.hist(ncbin[:-1], ncbin, weights = nchist / IC.widths, label=r"$\nu NC$", color="brown", histtype="step")
    ax1.hist(ncbarbin[:-1], ncbarbin, weights = ncbarhist / IC.widths, label=r"$\overline{\nu} NC$", color="brown",\
                            linestyle = "--", histtype="step")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
    ax1.set_xlim(IC.lo, IC.hi)
    ax1.set_ylabel("Effective Volume [m^3]")
    ax1.grid(True)
    ax1.legend(loc = 2)
    ax1.title.set_text("IceCube Upgrade Deepcore")

    # then plot ORCA
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.evol , label=r"$\nu_eCC$", color="red", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.ebarvol , label=r"$\overline{\nu}_eCC$", \
                            linestyle = '--', color="red", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.muvol , label=r"$\nu_{\mu}CC$", color="blue", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.mubarvol , label=r"$\overline{\nu}_{\mu}CC$", \
                            linestyle = '--', color="blue", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.tauvol , label=r"$\nu_{\tau}CC$", color="green", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.taubarvol , label=r"$\overline{\nu}_{\tau}CC$",\
                            color="green", linestyle = "--", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.ncvol , label=r"$\nu NC$", color="brown", histtype="step")
    ax2.hist(ORCA.e_bin[:-1], ORCA.e_bin, weights = ORCA.ncbarvol , label=r"$\overline{\nu} NC$", color="brown",\
                            linestyle = "--", histtype="step")
    ax2.set_xscale("log")
    # ax.set_yscale("log")
    ax2.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
    ax2.set_ylabel("Effective Volume [m^3]")
    ax2.set_xlim(1, 50)
    ax2.grid(True)
    ax2.legend(loc = 2)
    ax2.title.set_text("ORCA")

    fig.suptitle("IceCube Upgrade DeepCore and ORCA Effective Volumes")
    # plt.show()
    plt.savefig("./paper_plots/Effective_Volumes.png")
    plt.close()

effective_volumes()