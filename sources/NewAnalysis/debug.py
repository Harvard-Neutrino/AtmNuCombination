import numpy as np 
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import time
from matplotlib import pyplot as plt

from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext

import util
from params import *
import classdef as cl 

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

def IC_num_events():
    # first check number of total events for 3 years for IC
    sim = cl.Simulation(pd.read_csv(input_file), ORCA = False)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0
    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    num += analysis.bf_weights[i][j][k][l]

    print("IceCube upgrade predicted events per year is ", num / 3)
    return

def ORCA_num_events():
    # first check number of total events for 3 years for IC
    sim = cl.Simulation(pd.read_csv(input_file), ORCA = True)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0
    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    num += analysis.bf_weights[i][j][k][l]

    print("ORCA upgrade predicted events per year is ", num / 3)
    return

# this is to check what happened with the reweighting

# first check number of events

# first the control group
def control():
    sim = cl.Simulation(pd.read_csv("ORCA_control.csv"), ORCA = False)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0

    # also setup a overall propagated weight
    rated_weight = np.zeros_like(sim.W_mc)

    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    rated_weight[l] += analysis.bf_weights[i][j][k][l]
                    if np.cos(analysis.simulation.C_tr[l]) < 0:
                        num += analysis.bf_weights[i][j][k][l]

    print("ORCA upgrade predicted events per year is ", num / 3)
    print(rated_weight[:20])
    # now produce a histogram based on the rated weights and distributions
    track_mask = sim.pid == 1
    cas_mask = sim.pid == 0
    fig, [ax, ax2] = plt.subplots(1, 2, constrained_layout=True, figsize = (25, 12))
    x = np.logspace(np.log10(1.85), np.log10(54), 101) # this is energy
    y = np.linspace(-1, 0, 101) # this is zenith
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(sim.E_tr[track_mask], np.cos(sim.C_tr[track_mask]), \
                    bins=(x, y), weights = rated_weight[track_mask])
    Z2, xedges, yedges = np.histogram2d(sim.E_tr[cas_mask], np.cos(sim.C_tr[cas_mask]), \
                    bins=(x, y), weights = rated_weight[cas_mask])
    # attempt to manually normalize column

    im = ax.pcolor(X, Y, Z.T)
    im2 = ax2.pcolor(X, Y, Z2.T)

    ax.set_xlim(1.85, 54)
    ax.set_ylim(-1, 0)
    ax2.set_xlim(1.85, 54)
    ax2.set_ylim(-1, 0)
    plt.suptitle("Event Distribution Under No Reweighting")
    ax.title.set_text('Tracks')
    ax2.title.set_text('Cascades')
    ax.set_xlabel("True Energy")
    ax2.set_xlabel("True Energy")
    ax.set_ylabel("Cos True Zenith")
    fig.colorbar(im, orientation = "vertical")
    ax.set_xscale("log")
    ax2.set_xscale("log")
    # ax.set_yscale("log")
    # plt.legend()
    # plt.show()
    plt.savefig("./Plots/Event_Distribution_No_Reweighting")
    plt.close()

def control_reco():
    sim = cl.Simulation(pd.read_csv("ORCA_control.csv"), ORCA = False)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0

    # also setup a overall propagated weight
    rated_weight = np.zeros_like(sim.W_mc)

    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    rated_weight[l] += analysis.bf_weights[i][j][k][l]
                    if np.cos(analysis.simulation.C_tr[l]) < 0:
                        num += analysis.bf_weights[i][j][k][l]

    print("ORCA upgrade predicted events per year is ", num / 3)
    print(rated_weight[:20])
    # now produce a histogram based on the rated weights and distributions
    track_mask = sim.pid == 1
    cas_mask = sim.pid == 0
    fig, [ax, ax2] = plt.subplots(1, 2, constrained_layout=True, figsize = (25, 12))
    x = np.logspace(np.log10(1.85), np.log10(54), 101) # this is energy
    y = np.linspace(-1, 0, 101) # this is zenith
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(sim.E_re[track_mask], np.cos(sim.C_tr[track_mask]), \
                    bins=(x, y), weights = rated_weight[track_mask])
    Z2, xedges, yedges = np.histogram2d(sim.E_re[cas_mask], np.cos(sim.C_tr[cas_mask]), \
                    bins=(x, y), weights = rated_weight[cas_mask])
    # attempt to manually normalize column

    im = ax.pcolor(X, Y, Z.T)
    im2 = ax2.pcolor(X, Y, Z2.T)

    ax.set_xlim(1.85, 54)
    ax.set_ylim(-1, 0)
    ax2.set_xlim(1.85, 54)
    ax2.set_ylim(-1, 0)
    plt.suptitle("Event Distribution (Reco) Under No Reweighting")
    ax.title.set_text('Tracks')
    ax2.title.set_text('Cascades')
    ax.set_xlabel("True Energy")
    ax2.set_xlabel("True Energy")
    ax.set_ylabel("Cos True Zenith")
    fig.colorbar(im, orientation = "vertical")
    ax.set_xscale("log")
    ax2.set_xscale("log")
    # ax.set_yscale("log")
    # plt.legend()
    # plt.show()
    plt.savefig("./Plots/Event_Distribution_reco_No_Reweighting")
    plt.close()

# then the reweighting
def reweight():
    sim = cl.Simulation(pd.read_csv("ORCA_only_reweight.csv"), ORCA = False)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0

    # also setup a overall propagated weight
    rated_weight = np.zeros_like(sim.W_mc)

    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    rated_weight[l] += analysis.bf_weights[i][j][k][l]
                    if np.cos(analysis.simulation.C_tr[l]) < 0:
                        num += analysis.bf_weights[i][j][k][l]

    print("ORCA upgrade predicted events per year is ", num / 3)
    print(rated_weight[:20])
    # now produce a histogram based on the rated weights and distributions
    track_mask = sim.pid == 1
    cas_mask = sim.pid == 0
    fig, [ax, ax2] = plt.subplots(1, 2, constrained_layout=True, figsize = (25, 12))
    x = np.logspace(np.log10(1.85), np.log10(54), 101) # this is energy
    y = np.linspace(-1, 0, 101) # this is zenith
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(sim.E_tr[track_mask], np.cos(sim.C_tr[track_mask]), \
                    bins=(x, y), weights = rated_weight[track_mask])
    Z2, xedges, yedges = np.histogram2d(sim.E_tr[cas_mask], np.cos(sim.C_tr[cas_mask]), \
                    bins=(x, y), weights = rated_weight[cas_mask])
    # attempt to manually normalize column

    im = ax.pcolor(X, Y, Z.T)
    im2 = ax2.pcolor(X, Y, Z2.T)

    ax.set_xlim(1.85, 54)
    ax.set_ylim(-1, 0)
    ax2.set_xlim(1.85, 54)
    ax2.set_ylim(-1, 0)
    plt.suptitle("Event Distribution with Reweighting")
    ax.title.set_text('Tracks')
    ax2.title.set_text('Cascades')
    ax.set_xlabel("True Energy")
    ax2.set_xlabel("True Energy")
    ax.set_ylabel("Cos True Zenith")
    fig.colorbar(im, orientation = "vertical")
    ax.set_xscale("log")
    ax2.set_xscale("log")
    # ax.set_yscale("log")
    # plt.legend()
    # plt.show()
    plt.savefig("./Plots/Event_Distribution_with_Reweighting")
    plt.close()

def reweight_reco():
    sim = cl.Simulation(pd.read_csv("ORCA_only_reweight.csv"), ORCA = False)
    bf_fluxes = util.bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp)
    fluxes = bf_fluxes # dont care about chi squared for now
    analysis = cl.Analysis(sim, bf_fluxes, fluxes)
    util.get_all_weights(analysis, cl.PointType.BestFit)
    # util.get_all_weights(analysis, cl.PointType.Physical)
    # now some up all the propagated weights
    num = 0

    # also setup a overall propagated weight
    rated_weight = np.zeros_like(sim.W_mc)

    for i in range(2):
        for j in range(2):
            for k in range(3):
                for l in range(analysis.bf_weights.shape[-1]):
                    rated_weight[l] += analysis.bf_weights[i][j][k][l]
                    if np.cos(analysis.simulation.C_tr[l]) < 0:
                        num += analysis.bf_weights[i][j][k][l]

    print("ORCA upgrade predicted events per year is ", num / 3)
    print(rated_weight[:20])
    # now produce a histogram based on the rated weights and distributions
    track_mask = sim.pid == 1
    cas_mask = sim.pid == 0
    fig, [ax, ax2] = plt.subplots(1, 2, constrained_layout=True, figsize = (25, 12))
    x = np.logspace(np.log10(1.85), np.log10(54), 101) # this is energy
    y = np.linspace(-1, 0, 101) # this is zenith
    X, Y = np.meshgrid(x, y)
    Z, xedges, yedges = np.histogram2d(sim.E_re[track_mask], np.cos(sim.C_tr[track_mask]), \
                    bins=(x, y), weights = rated_weight[track_mask])
    Z2, xedges, yedges = np.histogram2d(sim.E_re[cas_mask], np.cos(sim.C_tr[cas_mask]), \
                    bins=(x, y), weights = rated_weight[cas_mask])
    # attempt to manually normalize column

    im = ax.pcolor(X, Y, Z.T)
    im2 = ax2.pcolor(X, Y, Z2.T)

    ax.set_xlim(1.85, 54)
    ax.set_ylim(-1, 0)
    ax2.set_xlim(1.85, 54)
    ax2.set_ylim(-1, 0)
    plt.suptitle("Event Distribution (Reco) with Reweighting")
    ax.title.set_text('Tracks')
    ax2.title.set_text('Cascades')
    ax.set_xlabel("True Energy")
    ax2.set_xlabel("True Energy")
    ax.set_ylabel("Cos True Zenith")
    fig.colorbar(im, orientation = "vertical")
    ax.set_xscale("log")
    ax2.set_xscale("log")
    # ax.set_yscale("log")
    # plt.legend()
    # plt.show()
    plt.savefig("./Plots/Event_Distribution_reco_with_Reweighting")
    # plt.close()
# control()
# reweight()
reweight_reco()
# control_reco()


