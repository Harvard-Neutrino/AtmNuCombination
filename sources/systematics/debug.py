import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import nuSQuIDS as nsq
import nuflux
import time
# import seaborn as sns

from params import *
import chisq as chi
import util
import plotting
import propagate as prop


# start = time.time()
# truth0 = chi.get_truth(0)
# end = time.time()
# print("truth took {}".format(end - start))
# # truth1 = chi.get_truth(1)
# # print(truth0)
# start = time.time()
# rate_weight = prop.propagate(theta23, m31)
# end = time.time()

# print("propagate took {}".format(end - start))
# # print(rate_weight[0:3])

# chisqval = chi.min_chisq(rate_weight, truth0, 0) #+ chi.min_chisq(rate_weight, truth1, 1)

# plotting.plot_contour("1108_sys_2D_trial_2")
plotting.plot_profile(1, "1109_m_profile_norm_vs_no_sys")


# no_syst = util.read_output(dir_name = "../NoSystematics/1109_no_sys_the_profile")[0]
# profile = util.read_output(dir_name = "./1109_theta_profile_norm_min/")[0]

# for i in range(len(t23l)):
#     print(np.sin(t23l[i])**2, no_syst[i], profile[i])
