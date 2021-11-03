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

# print(util.read_output())
# print(len(t23l), t23l)
# print(len(m31l), m31l)
plotting.plot_contour("chisq_contour_1102TrialRun3")