import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux
import time
pd.options.mode.chained_assignment = None  # default='warn'

from params import *
import propagate as prop

def add_systematics(W_r, delta, gamma, top):
    input_data["rate_weight"] = W_r #temp[:]

    # apply nu/nubar ratio and slope
    for i in range(len(input_data["rate_weight"])):
        # print("1")
        if input_data["pdg"][i] > 0:
            # print("2")
            input_data["rate_weight"][i] *= delta
        if input_data["pdg"][i] < 0:
            # print("3")
            input_data["rate_weight"][i] *= (2 - delta)
        # print(input_data["true_energy"])
        temp1 = input_data["true_energy"][i] / E_0
        # print(temp1)
        # print("4")
        temp2 = temp1**gamma
        # print(temp2)
        # print("5")
        input_data['rate_weight'][i] *= temp2
        

            
    if top == 0:
        hist_truth, _, _ = np.histogram2d(x = input_data["reco_energy"][cascade_mask], y = np.cos(input_data["reco_zenith"][cascade_mask]), bins = [energy_bins_fine, cos_bin_plot], weights = input_data["rate_weight"][cascade_mask])
    elif top == 1:
        hist_truth, _, _ = np.histogram2d(x = input_data["reco_energy"][track_mask], y = np.cos(input_data["reco_zenith"][track_mask]), bins = [energy_bins_fine, cos_bin_plot], weights = input_data["rate_weight"][track_mask])
    else:
        hist_truth, _, _ = np.histogram2d(x = input_data["reco_energy"], y = np.cos(input_data["reco_zenith"]), bins = [energy_bins_fine, cos_bin_plot], weights = input_data["rate_weight"])


    return hist_truth