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

'''
# This part is testing the time used 
start = time.time()
truth0 = chi.get_truth(0)
end = time.time()
print("truth took {}".format(end - start))
# truth1 = chi.get_truth(1)
# print(truth0)
start = time.time()
rate_weight = prop.propagate(theta23, m31)
end = time.time()

print("propagate took {}".format(end - start))
# print(rate_weight[0:3])

chisqval = chi.min_chisq(rate_weight, truth0, 0) #+ chi.min_chisq(rate_weight, truth1, 1)
'''

print(chi.get_chisq(theta23, m31, chi.get_truth(0), 0))