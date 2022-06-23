import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit
import math
# from Effective import ICEffectiveAnalysis as ICEff
# from Effective import ORCAEffectiveAnalysis as ORCAEff
import digitalizer as dgt
import util
from util import gaussian
# import analyze as anl
from params import *

input_file = pd.read_csv("0622ORCA_no_morph_with_LE.csv")
original = pd.read_csv("neutrino_mc.csv")
pid = input_file["pid"]
e_true = input_file["true_energy"]
e_reco = input_file["reco_energy"]
zen_true = input_file["true_zenith"]
zen_reco = input_file["reco_zenith"]
weights = input_file["weight"]

# cnt0 = 0
# cnt1 = 0
# cnt2 = 0
# for i in range(len(weights)):
#     if pid[i] == 2 and weights[i] != 0:
#         print(i)
#     if pid[i] == 0 and weights[i] == 0:
#         print(i)
#     if pid[i] == 0:
#         cnt0 += 1
#     if pid[i] == 1:
#         cnt1 += 1
#     if pid[i] == 2:
#         cnt2 += 1
# print(cnt0, cnt1, cnt2)


# mask = []
# for i in range(len(original["weight"])):
# 	if original["true_energy"][i] >= 54 or original["true_energy"][i] <= 1:
# 		mask.append(False)
# 	else:
# 		mask.append(True)

# input_MC = original
# df_array = np.ndarray((17, len(input_MC["pdg"][mask])))
# df_array[0] = input_MC["pdg"][mask]
# df_array[1] = input_MC["pid"][mask]
# df_array[2] = input_MC["true_energy"][mask]
# df_array[3] = input_MC["true_zenith"][mask]
# df_array[4] = input_MC["true_azimuth"][mask]
# df_array[5] = input_MC["reco_energy"][mask]
# df_array[6] = input_MC["reco_zenith"][mask]
# df_array[7] = input_MC["reco_azimuth"][mask]
# df_array[8] = input_MC["interaction_type"][mask]
# df_array[9] = input_MC["current_type"][mask]
# df_array[10] = input_MC["weight"][mask]
# df_array[11] = input_MC["xsec"][mask]
# df_array[12] = input_MC["dxsec"][mask]
# df_array[13] = input_MC["x"][mask]
# df_array[14] = input_MC["y"][mask]
# df_array[15] = input_MC["W"][mask]
# df_array[16] = input_MC["Q2"][mask]
# df = pd.DataFrame(df_array.T, columns = ["pdg", "pid", "true_energy", "true_zenith", \
# 	"true_azimuth", "reco_energy", "reco_zenith", "reco_azimuth", "interaction_type",\
# 	"current_type", "weight", "xsec", "dxsec", "x", "y", "W", "Q2"])


# df.to_csv("selected_ICMC")

