import sys
import os
import numpy as np

import chisq as chi
import util
import propagate as prop
import fisher
from params import *

# this file calculates a chi-squared value for a pair of physical parameters
# and then saves the physical parameters as well as the chi-sq into a file

idx = int(sys.argv[1])

# get the physical parameters
f1val, f2val = util.fisher_id_to_param(idx, fisher1, fisher2)

print("fisher1 and fisher2 are: {}, {}".format(f1val, f2val))

file_name = "{}.txt".format(idx)
complete_name = os.path.join(dir_name, file_name)

print("begin truth calculations")
truth0 = chi.get_truth(0)
truth1 = chi.get_truth(1)

if not fisher1_is_theta and not fisher2_is_theta:
    t23val = theta23
else:
    if fisher1_is_theta:
        t23val = f1val
    else:
        t23val = f2val

if not fisher1_is_m and not fisher2_is_m:
    m31val = m31
else:
    if fisher1_is_m:
        m31val = f1val
    else:
        m31val = f2val

print("begin calculating d2chi")
d2chival = fisher.d2chi_dtheta2(f1val, truth0, truth1)

saveres = np.array([idx, d2chiqval])

np.savetxt(complete_name, saveres)