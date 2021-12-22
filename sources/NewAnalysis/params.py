import numpy as np
import pandas as pd
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux
import sys
from scipy.optimize import minimize

units = nsq.Const()
Time = 3.0*365*24*60*60
meter_to_cm_sq = 1e4
Radian = np.pi / 180.

input_file = "neutrino_mc.csv"

interactions = False
neutrino_flavors = 3

E_min = 1.0*units.GeV
E_max = 1.0e3*units.GeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

theta12 = np.arcsin(np.sqrt(0.304))
theta13 = np.arcsin(np.sqrt(0.02221))
theta23 = np.arcsin(np.sqrt(0.570))
dm21 = 7.42e-5
dm31 = 2.517e-3
dcp = 0

Erec_min = 1
Erec_max = 10000