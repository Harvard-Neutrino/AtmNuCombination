import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux
import time

from params import *
import propagate as prop
import systematics as syst

def get_truth(top):
    W_r_bf = prop.propagate(theta23, m31)
    res, _ = syst.add_systematics(W_r_bf, top)
    return res

# Get chisq for the contour plot
def get_chisq(W_r, H, N0, H_bf, top):
    # Get the rated weight and energy bins for the given t23, m31 and truth

    # energy_bins, _ = syst.add_systematics(W_r, top)
    # energy_bins = H[:]

    chisq = 0
    for i in range(len(H)):
        # chisqplus = (N0 * energy_bins[i] - H_bf[i]) ** 2 /  H_bf[i]
        chisqplus = (N0 * H[i] - H_bf[i]) ** 2 /  H_bf[i]
        chisq += chisqplus
    # punishment for normalization
    normloss = ((N0 - 1) / SigmaN0) ** 2
    chisq += normloss
    return chisq

# Find the minimum chisq with systematics
def min_chisq(W_r, truth, top):

    # this is for testing if code works on cluster
    Hist, _ = syst.add_systematics(W_r, top)

    ls = []
    for n in N0l:
        ls.append(get_chisq(W_r, Hist, n, truth, top))
    minidx = ls.index(min(ls))
    minchi = ls[minidx]
    print(ls)
    return minchi