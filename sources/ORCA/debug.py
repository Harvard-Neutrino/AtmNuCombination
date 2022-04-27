import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit
import math

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

for i in range(len(pid)):
    if input_file["weight"][i] != 0:
        if zen_reco[i] == -1:
            print(i, " zen")
        if e_reco[i] == -1:
            print(i, " e")