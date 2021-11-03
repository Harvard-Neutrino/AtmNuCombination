import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux

from params import *
import propagate as prop

def add_systematics(theta23in, m31in, N0, top):
    rate_weight = prop.propagate(theta23in, m31in)
    # apply normalization
    if N0 != 1:
        for i in range(len(rate_weight)):
            rate_weight[i] *= N0
    
    input_data["rate_weight"] = rate_weight
    if top == 0:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"][cascade_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][cascade_mask])
    elif top == 1:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"][track_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][track_mask])
    else:
        energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"], bins = energy_bins_fine, weights = input_data["rate_weight"])
    
    return rate_weight, energy_hist_truth, energy_bins_truth