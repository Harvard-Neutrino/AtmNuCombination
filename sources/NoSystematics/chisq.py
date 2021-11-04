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

def get_truth(top):
    _, H_bf, _ = prop.propagate(theta23, m31)
    return res

# Get chisq for the contour plot
def get_chisq(t23, m31, H_bf, top):
    # Get the rated weight and energy bins for the given t23, m31 and truth
    _, H, _ = prop.propagate(t23, m31, top)
    chisq = 0
    for i in range(len(energy_bins)):
        chisqplus = (H[i] - H_bf[i]) ** 2 /  H_bf[i]
        chisq += chisqplus

    return chisq

