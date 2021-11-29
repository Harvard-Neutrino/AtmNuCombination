import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQuIDS as nsq
import nuflux
import seaborn as sns



# Set up mixing parameters as in nu-fit5.0
theta12 = np.arcsin(np.sqrt(0.304))
theta13 = np.arcsin(np.sqrt(0.02221))
theta23 = np.arcsin(np.sqrt(0.570))
m21 = 7.42e-5
m31 = 2.517e-3 # real Nu-Fit value
# m31 = 2.8e-3  # try some larger parameters

dir_name = "./1128_fisher_theta/"

# thetalist = True
# mlist = True

# Set theta23 numeric values to probe sensitivity
t23min = np.arcsin(np.sqrt(0.33))
t23max = np.arcsin(np.sqrt(0.67))
t23step = 0.001 * np.pi
# if thetalist:
t23l = np.arange(t23min, t23max + t23step, t23step)
# else:
#     t23l = np.array([theta23])
# t23l = np.array([0.44872923, theta23])

# if mlist:
# m31min = 2.20e-3
# m31max = 2.80e-3
# m31step = 0.01e-3
# else:
m31min = m31
m31max = m31
m31step = 1
m31l = np.arange(m31min, m31max + m31step, m31step)

# Set the fisher information plotting points lists
fisher_null = np.array([0])
fisher_t23l = np.arange(np.arcsin(np.sqrt(0.33)), np.arcsin(np.sqrt(0.67)) + 0.001 * np.pi, 0.001 * np.pi)
fisher_m31l = np.arange(2.20e-3, 2.80e-3 + 0.01e-3, 0.01e-3)

fisher1 = fisher_t23l
fisher2 = fisher_null
fisher1_is_theta = (fisher1.all() == fisher_t23l.all())
fisher1_is_m = (fisher1.all() == fisher_m31l.all())
fisher2_is_theta = (fisher2.all() == fisher_t23l.all())
fisher2_is_m = (fisher2.all() == fisher_m31l.all())


# set up systematics constants
E_0 = 20

# set up normalization range
SigmaN0 = 0.4
N0min = 0.9
N0max = 1.1
N0step = 0.001
N0l = np.arange(N0min, N0max + N0step, N0step)
# N0l = [1, 0.6, 1, 1.4]

# Define path to file (you may need to change this to match your system)
input_file = "neutrino_mc.csv"

# Load the file using pandas
input_data = pd.read_csv(input_file)

# Define masks to identify different neutrino flavors
nue_mask = (np.abs(input_data["pdg"]) == 12)
numu_mask = (np.abs(input_data["pdg"]) == 14)
nutau_mask = (np.abs(input_data["pdg"]) == 16)

# Define masks to identify interaction species
cascade_mask = (np.abs(input_data["pid"]) == 0)
track_mask = (np.abs(input_data["pid"]) == 1)

# Define masks to identify different flavor/interaction combinations.
nc_mask = input_data["current_type"] == 0
cc_mask = input_data["current_type"] == 1
nue_cc_mask = nue_mask & cc_mask
numu_cc_mask = numu_mask & cc_mask
nutau_cc_mask = nutau_mask & cc_mask
nue_nc_mask = nue_mask & nc_mask
numu_nc_mask = numu_mask & nc_mask
nutau_nc_mask = nutau_mask & nc_mask 
nue_cascade_mask = nue_mask & cascade_mask
nue_track_mask = nue_mask & track_mask
nutau_cascade_mask = nutau_mask & cascade_mask
nutau_track_mask = nutau_mask & track_mask
numu_cascade_mask = numu_mask & cascade_mask
numu_track_mask = numu_mask & track_mask

# Define some energy bins (used throughout this notebook)
energy_bins_fine = np.logspace(0., 3., num=31)
energy_bins_course = np.logspace(0., 2., num=11)

# Define the units and interaction mode
units = nsq.Const()
interactions = False

# Set propagation bins
E_min = 1*units.GeV
E_max = 1.0e3*units.GeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3

# Set the chi squared plotting bins limits
E_bin_min = 0
E_bin_max = 3
E_n_bins = 31
E_bin_plot = np.logspace(E_bin_min, E_bin_max, E_n_bins)

cos_bin_min = -1
cos_bin_max = 1
cos_n_bins = 11
cos_bin_plot = nsq.linspace(cos_bin_min, cos_bin_max, cos_n_bins)

theta_bin_min = np.pi
theta_bin_max = 2 * np.pi
theta_n_bins = 11
theta_bin_plot = nsq.linspace(theta_bin_min, theta_bin_max, theta_n_bins)


# Set up chi squared bins
bins = np.zeros((t23l.shape[0], E_n_bins))

