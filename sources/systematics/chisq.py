import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux

from params import *
import propagate as prop
import systematics as syst

def get_truth(top):
    rate_weight = prop.propagate(theta23, m31)
    res, _ = syst.add_systematics(rate_weight, 1, top)
    return res

# Get chisq for the contour plot
def get_chisq(W_r, N0, truth, top):
    # Get the rated weight and energy bins for the given t23, m31 and truth
    # rate_weight = prop.propagate(t23, m31)
    temprate2 = W_r[:]
    # print("get_chi before syst", W_r[0:3])
    energy_bins, _ = syst.add_systematics(temprate2, N0, top)
    # print("get_chi after syst", W_r[0:3])
    # print(N0, energy_bins)
    # rate_weight_truth, energy_hist_truth, energy_bins_truth = get_rated_weight_truth(top)
    chisq = 0
    for i in range(len(energy_bins)):
        # chisqplus = (energy_bins[i] - energy_hist_truth[i]) ** 2 /  energy_hist_truth[i]
        chisqplus = (energy_bins[i] - truth[i]) ** 2 /  truth[i]
        chisq += chisqplus
    # punishment for normalization
    normloss = ((N0 - 1) / SigmaN0) ** 2
    chisq += normloss
    return chisq

# Find the minimum chisq with systematics
def min_chisq(W_r, truth, top):
    ls = []
    for n in N0l:
        # print("min_chi loop", W_r[0:3])
        temprate = W_r[:]
        ls.append(get_chisq(temprate, n, truth, top))
    print(ls)
    minidx = ls.index(min(ls))
    minchi = ls[minidx]
    return minchi