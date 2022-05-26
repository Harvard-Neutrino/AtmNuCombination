import numpy as np 
import pandas as pd

input_scale = "./ORCA_Results/scale.png"
input_track = "./ORCA_Results/migrationTracks.png"
input_cascade = "./ORCA_Results/migrationCascades.png"

input_MC = pd.read_csv("neutrino_mc.csv")
savename = "ORCA_only_no_tau.csv"

x_bins = np.logspace(np.log10(1.85), np.log10(54), 23)

# turn on or off randomize morphology
rand_morph = True

# turn on or off pessimistic ereco for intermediate class events
pess = False

# ratio of track vs cascade reco of intermediate events
if not pess:
	pess_ereco = 2
elif pess:
	pess_ereco = 3

# turn on or off random zenith reco
rand_zen = True

# turn on or off random energy reco
rand_energy = True

# turn on or off MC weights
reweight = True

# turn on or off tau neutrinos reweighting
reweight_tau = False

# turn on or off two log norm fit instead of one log norm fit !!!!! CURRENTLY DOES NOT WORK !!!!!
two_gaus = False
