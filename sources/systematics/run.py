import sys
import os

import chisq as chi
import util
from params import *

# this file calculates a chi-squared value for a pair of physical parameters
# and then saves the physical parameters as well as the chi-sq into a file

idx = int(sys.argv[1])
dir_name = "./1102TrialRun/"

# get the physical parameters
t23val, m31val = util.id_to_param(idx)

file_name = "{}_{}.txt".format(t23val, m31val)
complete_name = os.path.join(dir_name, file_name)

truth0 = chi.get_truth(0)
truth1 = chi.get_truth(1)
chisqval = chi.get_chisq(t23val, m31val, truth0, 0) + chi.get_chisq(t23val, m31val, truth1, 1)

file1 = open(complete_name, "a")

L = ["{}\n".format(idx), "{}\n".format(t23val), "{}\n".format(m31val), str(chisqval)]
file1.writelines(L) 

file1.close()

