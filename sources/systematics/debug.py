import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux
# import seaborn as sns

from params import *
import chisq as chi
import util
import plotting
import propagate as prop



truth0 = chi.get_truth(0)
# truth1 = chi.get_truth(1)
# print(truth0)

rate_weight = prop.propagate(theta23, m31)
# print(rate_weight[0:3])

chisqval = chi.min_chisq(rate_weight, truth0, 0) #+ chi.min_chisq(rate_weight, truth1, 1)