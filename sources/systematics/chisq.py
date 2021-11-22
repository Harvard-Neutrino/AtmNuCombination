import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux
import time
import scipy as scp

from params import *
import propagate as prop
import systematics as syst

def get_truth(top):
    W_r_bf = prop.propagate(theta23, m31)
    res = syst.add_systematics(W_r_bf, 1, 1, top)
    return res

# Get chisq for the contour plot
def get_chisq(W_r, N0, delta, gamma, H_bf, top):
    # Get the rated weight and energy bins for the given t23, m31 and truth

    H = syst.add_systematics(W_r, delta, gamma, top)

    chisq = 0
    (rows, cols) = H.shape
    for i in range(rows):
        for j in range(cols):
            chisqplus = (N0 * H[i][j] - H_bf[i][j]) ** 2 /  H_bf[i][j]
            chisq += chisqplus
    # punishment for normalization
    normloss = ((N0 - 1) / SigmaN0) ** 2
    deltaloss = ((delta - 1) / SigmaN0) ** 2 # use sigma_delta = sigma_n0
    gammaloss = ((gamma - 1) / SigmaN0) ** 2 # use sigma_gamma = sigma_n0
    chisq += normloss
    return chisq

Nfeval = 1

# Find the minimum chisq with systematics
def min_chisq(W_r, truth, top):

    # Hist = syst.add_systematics(W_r, top)

    # define the function with only 1 variable to go through the minimization
    def to_min(syst):
        # syst := np.array([N0, delta, gamma])
        return get_chisq(W_r, syst[0], syst[1], syst[2], truth, top)
    
    # define a callback function to visualize the results
    def callbackF(syst):
        global Nfeval
        print('{0:4d}, {1: 3.6f}, {2: 3.6f}, {3: 3.6f}, {3: 3.6f}'.format\
                                (Nfeval, syst[0], syst[1], syst[2], to_min(syst)))
        Nfeval += 1
    
    res = scp.optimize.minimize(to_min, np.array([1, 1, 1]), callback = callbackF, options={'disp': True, 'maxiter': 30})
    # res = scp.optimize.minimize(to_min, np.array([1]), options={'disp': True})


    return res.fun