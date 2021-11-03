import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux
import time
pd.options.mode.chained_assignment = None  # default='warn'

import params as pm
import propagate as prop

def add_systematics(W_r, top):
    # apply normalization

    temp = W_r[:]
    pm.input_data["rate_weight"] = temp[:]
            
    if top == 0:
        energy_hist_truth, energy_bins_truth = np.histogram(pm.input_data["reco_energy"][pm.cascade_mask], bins = pm.energy_bins_fine, weights = pm.input_data["rate_weight"][pm.cascade_mask])
    elif top == 1:
        energy_hist_truth, energy_bins_truth = np.histogram(pm.input_data["reco_energy"][pm.track_mask], bins = pm.energy_bins_fine, weights = pm.input_data["rate_weight"][pm.track_mask])
    else:
        energy_hist_truth, energy_bins_truth = np.histogram(pm.input_data["reco_energy"], bins = pm.energy_bins_fine, weights = pm.input_data["rate_weight"])
    
    return energy_hist_truth, energy_bins_truth