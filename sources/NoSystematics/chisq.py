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
    _, H_bf = prop.propagate(theta23, m31, top)
    print(H_bf)
    return H_bf

# Get chisq for the contour plot
def get_chisq(t23, m31, H_bf, top):
    # Get the rated weight and energy bins for the given t23, m31 and truth
    _, H = prop.propagate(t23, m31, top)
    chisq = 0
    (rows, cols) = H.shape
    for i in range(rows):
        # print(i)
        for j in range(cols):
            # print(j)
            chisqplus = (H[i][j] - H_bf[i][j]) ** 2 /  H_bf[i][j]
            chisq += chisqplus

    return chisq

