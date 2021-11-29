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
import chisq as chi

def d2chi_dtheta2(theta, truth0, truth1):
    def chieval(t23val):
        rate_weight = prop.propagate(t23val, m31)

        chisqval = chi.min_chisq(rate_weight, truth0, 0) + chi.min_chisq(rate_weight, truth1, 1)
        return chisqval
    return scp.misc.derivative(chieval, theta, dx = 1e-7, n = 2, order = 7)