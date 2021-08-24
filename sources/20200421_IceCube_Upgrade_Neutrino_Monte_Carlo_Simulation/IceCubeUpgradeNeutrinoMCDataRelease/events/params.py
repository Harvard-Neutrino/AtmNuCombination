import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq
import nuflux
import seaborn as sns

# Set up mixing parameters as in nu-fit5.0
theta12 = np.arcsin(np.sqrt(0.304))
theta13 = np.arcsin(np.sqrt(0.02221))
theta23 = np.arcsin(np.sqrt(0.570))
m21 = 7.42e-5
m31 = 2.517e-3


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
energy_bins_fine = np.logspace(0., 2., num=21)
energy_bins_course = np.logspace(0., 2., num=11)

# Define the units and interaction mode
units = nsq.Const()
interactions = False

# Set propagation bins
E_min = 10.0*units.GeV
E_max = 1.0e3*units.GeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3

# Set theta23 numeric values to probe sensitivity
t23min = np.sqrt(np.arcsin(0.3))
t23max = np.sqrt(np.arcsin(0.8))
t23step = 0.0025 * np.pi
t23l = np.arange(t23min, t23max + t23step, t23step)
# t23l = np.array([0.44872923, theta23])

m31min = 2e-3
m31max = 4e-3
m31step = 0.1e-3
m31l = np.arange(m31min, m31max + m31step, m31step)

# Set the chi squared plotting bins limits
E_bin_min = 0
E_bin_max = 2
E_n_bins = 20
E_bin_plot = np.logspace(E_bin_min, E_bin_max, E_n_bins)

cos_bin_min = -1
cos_bin_max = 1
cos_n_bins = 10
cos_bin_plot = nsq.linspace(cos_bin_min, cos_bin_max, cos_n_bins)

theta_bin_min = np.pi
theta_bin_max = 2 * np.pi
theta_n_bins = 10
theta_bin_plot = nsq.linspace(theta_bin_min, theta_bin_max, theta_n_bins)


# Set up chi squared bins
bins = np.zeros((t23l.shape[0], E_n_bins))



# # Set up some more masks for plotting bias (resolution + unfolded)
# reso_mask0 = input_data["true_energy"] >= 10^(1.0) & input_data["true_energy"] < 10^(1.1)
# reso_mask1 = input_data["true_energy"] >= 10^(1.1) & input_data["true_energy"] < 10^(1.2)
# reso_mask2 = input_data["true_energy"] >= 10^(1.2) & input_data["true_energy"] < 10^(1.3)
# reso_mask3 = input_data["true_energy"] >= 10^(1.3) & input_data["true_energy"] < 10^(1.4)
# reso_mask4 = input_data["true_energy"] >= 10^(1.4) & input_data["true_energy"] < 10^(1.5)
# reso_mask5 = input_data["true_energy"] >= 10^(1.5) & input_data["true_energy"] < 10^(1.6)
# reso_mask6 = input_data["true_energy"] >= 10^(1.6) & input_data["true_energy"] < 10^(1.7)
# reso_mask7 = input_data["true_energy"] >= 10^(1.7) & input_data["true_energy"] < 10^(1.8)
# reso_mask8 = input_data["true_energy"] >= 10^(1.8) & input_data["true_energy"] < 10^(1.9)
# reso_mask9 = input_data["true_energy"] >= 10^(1.9) & input_data["true_energy"] < 10^(2.0)
