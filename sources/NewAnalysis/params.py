import numpy as np
import pandas as pd
import nuSQuIDS as nsq
import nuflux
import sys
from scipy.optimize import minimize

units = nsq.Const()
Time = 3.0*365*24*60*60
meter_to_cm_sq = 1e4
Radian = np.pi / 180.

######################################################
NT23 = 24
SqT23_min = 0.305
SqT23_max = 0.705
SqT23 = np.linspace(SqT23_min, SqT23_max, NT23+1, endpoint = True)
SqT23BF = SqT23[17]
T23 = np.arcsin(np.sqrt(SqT23))
T23BF = T23[16]

NDM31 = 24
DM31_max = 3.0e-3
DM31_min = 2.0e-3
DM31 = np.linspace(DM31_min, DM31_max, NDM31+1, endpoint = True)
DM31BF = DM31[12]

NDCP = 19
DCP_min = 0
DCP_max = 2*np.pi
DCP = np.linspace(DCP_min, DCP_max, NDCP+1, endpoint = True)
DCPBF = DCP[13]
######################################################

input_file = "neutrino_mc.csv"

# define some nusquids parameters
interactions = False
neutrino_flavors = 3

# define some propagation bins
E_min = 1.0*units.GeV
E_max = 1.0e3*units.GeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
num_cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max, num_cth_nodes)

# define some physical parameters
theta12 = np.arcsin(np.sqrt(0.304))
theta13 = np.arcsin(np.sqrt(0.02221))
theta23 = np.arcsin(np.sqrt(0.570))
dm21 = 7.42e-5
dm31 = 2.517e-3
dcp = 4.299


# define some reconstruction bins
Erec_min = 1
Erec_max = 10000
NErec = 40
erec = np.logspace(np.log10(Erec_min), np.log10(Erec_max), NErec+1, endpoint = True)
dlter = (np.log10(Erec_max) - np.log10(Erec_min))/NErec

Crec_min = -1
Crec_max = 1
Ncrec = 10
crec = np.linspace(Crec_min, Crec_max, Ncrec+1, endpoint = True)
dltcr = (Crec_max - Crec_min)/Ncrec

# define some systematics parameters
# Energy slope
E0 = 10
gamma_bf = 0
sig_gamma = 0.2

# normalization
N_bf = 1
sig_N = 0.4

# nu-nubar ratio
delta_bf = 1
sig_delta = 0.02

# horizontal-vertical ratio
hv_bf = 0
sig_hv = 0.2

# electron-muon ratio
eps_bf = 1
sig_eps = 0.05


