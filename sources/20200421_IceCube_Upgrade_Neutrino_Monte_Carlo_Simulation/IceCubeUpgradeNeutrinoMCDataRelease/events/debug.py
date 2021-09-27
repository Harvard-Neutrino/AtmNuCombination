import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQuIDS as nsq
import nuflux
import seaborn as sns

import sensitivity as sst
from params import *

top = 0
rate_weight_truth, energy_hist_truth, energy_bins_truth = sst.get_rated_weight_truth(top)

def circle(x, y):
    return x ** 2 + y ** 2

def contour():
    fig, ax = plt.subplots()
    x = np.arange(-5, 5, 0.1)
    y = np.arange(-5, 5, 0.1)
    X, Y = np.meshgrid(x, y)

    Z = circle(X, Y)
    print(Z)

    Z2 = np.ndarray((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            Z2[i][j] = circle(x[i], x[j])
    print(Z2)

    axim = ax.contour(X,Y,Z2,levels=[1e0,1e1,1e2], cmap=plt.cm.jet)
    cb   = fig.colorbar(axim)
    plt.show()

# sst.get_chisq(theta23, m31, energy_hist_truth, top)
# t23l = [theta23]
sst.get_t23_chi_profile()

