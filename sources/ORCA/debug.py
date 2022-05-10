import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit
import math
from Effective import ICEffectiveAnalysis as ICEff
from Effective import ORCAEffectiveAnalysis as ORCAEff
import digitalizer as dgt
import util
from util import gaussian
# import analyze as anl
from params import *

input_file = pd.read_csv("ORCA.csv")
original = pd.read_csv("neutrino_mc.csv")
pid = input_file["pid"]
e_true = input_file["true_energy"]
e_reco = input_file["reco_energy"]
zen_true = input_file["true_zenith"]
zen_reco = input_file["reco_zenith"]
weights = input_file["weight"]

# ICVol = ICEff(1, 50, 51)
# ICVol.computeArea()
# ICVol.computeVolume()
# ORCAVol = ORCAEff()
# ORCAVol.computeArea()
# ORCAVol.computeVolume()
# ICe, ICeb, ICmu, ICmub, ICtau, ICtaub, ICnc, ICncb = ICVol.returnVol()
# ORe, OReb, ORmu, ORmub, ORtau, ORtaub, ORnc, ORncb = ORCAVol.returnVol()

# e = ORe/ICe
# mu = ORmu/ICmu
# tau = ORtau/ICtau
# nc = ORnc/ICnc
# eb = OReb/ICeb
# mub = ORmub/ICmub
# taub = ORtau/ICtaub
# ncb = ORnc/ICncb

# for i in range(17):
#     tau[i] = 0
#     taub[i] = 0

# for i in range(19):
#     if tau[i] > 30:
#         tau[i] = 30
#     if taub[i] > 30:
#         taub[i] = 30

# print(e)
# print(mu)
# print(tau)
# print(nc)
# print(eb)
# print(mub)
# print(taub)
# print(ncb)

# sorted = input_file.sort_values('weight', ascending = False)
# print(sorted['weight'][:10])
# print(sorted['true_energy'][:10])

# sorted2 = original.sort_values('weight', ascending = False)
# print(sorted2['weight'][:10])
# print(sorted2['true_energy'][:10])

# for i in range(len(pid)):
#     if input_file["weight"][i] != 0:
#         if zen_reco[i] == -1:
#             print(i, " zen")
#         if e_reco[i] == -1:
#             print(i, " e")