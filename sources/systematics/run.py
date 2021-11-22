import sys
import os
import numpy as np

import chisq as chi
import util
import propagate as prop
from params import *

# this file calculates a chi-squared value for a pair of physical parameters
# and then saves the physical parameters as well as the chi-sq into a file

idx = int(sys.argv[1])

# get the physical parameters
t23val, m31val = util.id_to_param(idx)

print("t23 and m31 are: {}, {}".format(t23val, m31val))

file_name = "{}.txt".format(idx)
complete_name = os.path.join(dir_name, file_name)

truth0 = chi.get_truth(0)
truth1 = chi.get_truth(1)

rate_weight = prop.propagate(t23val, m31val)

print("begin the min_chisq")

chisqval = chi.min_chisq(rate_weight, truth0, 0) + chi.min_chisq(rate_weight, truth1, 1)

# print("the final chi square is: ", chisqval)
# exit(0)

saveres = np.array([idx, chisqval])

np.savetxt(complete_name, saveres)

# file1 = open(complete_name, "a")

# L = ["{}\n".format(idx), "{}\n".format(t23val), "{}\n".format(m31val), str(chisqval)]
# file1.writelines(L) 

# file1.close()

