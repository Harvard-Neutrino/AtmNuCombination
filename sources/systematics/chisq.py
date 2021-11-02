import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux

from params import *
import propagate as prop

def get_truth(top):
    _, res, _ = prop.propagate(theta23, m31, top)
    return res

# Get chisq for the contour plot
def get_chisq(t23, m31, truth, top):
    # Get the rated weight and energy bins for the given t23, m31 and truth
    _, energy_bins, _ = prop.propagate(t23, m31, top)
    # rate_weight_truth, energy_hist_truth, energy_bins_truth = get_rated_weight_truth(top)
    chisq = 0
    for i in range(len(energy_bins)):
        # chisqplus = (energy_bins[i] - energy_hist_truth[i]) ** 2 /  energy_hist_truth[i]
        chisqplus = (energy_bins[i] - truth[i]) ** 2 /  truth[i]
        chisq += chisqplus
    return chisq
