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
import systematics as syst

# util.read_output(dir_name = "./1128_theta_3_sys_min/")[0]
t23val = 0.7910104958884872
m23val = 0.002517

truth0 = chi.get_truth(0)
truth1 = chi.get_truth(1)

W_r = prop.propagate(t23val, m31val)

print(chi.get_chisq(W_r, N0, delta, gamma, truth0, 0))

'''

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

# plotting.plot_contour("1109_countour_norm_min")
# plotting.plot_profile(0, "1117_theta_norm+delta_min")

truth0 = chi.get_truth(0)
truth1 = chi.get_truth(1)

W_r = prop.propagate(np.arcsin(np.sqrt(0.5)), m31)

print(chi.get_chisq(W_r, 1, 1, truth0, 0) + chi.get_chisq(W_r, 1, 1, truth1, 1))

# no_syst = util.read_output(dir_name = "../NoSystematics/1109_no_sys_the_profile")[0]
# profile = util.read_output(dir_name = "./1109_theta_profile_norm_min/")[0]

# for i in range(len(t23l)):
#     print(np.sin(t23l[i])**2, no_syst[i], profile[i])
'''