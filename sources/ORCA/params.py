import numpy as np 
import pandas as pd

input_scale = "./ORCA_Results/scale.png"
input_track = "./ORCA_Results/migrationTracks.png"
input_cascade = "./ORCA_Results/migrationCascades.png"

input_MC = pd.read_csv("neutrino_mc.csv")

x_bins = np.logspace(np.log10(1.85), np.log10(54), 23)

# turn on or off randomize morphology
rand_morph = True

# turn on or off pessimistic ereco for intermediate class events
pess = False

if not pess:
	pess_ereco = 2
elif pess:
	pess_ereco = 3

# turn on or off random zenith reco
rand_zen = False

# turn on or off random energy reco
rand_energy = False

# turn on or off MC weights
reweight = False

# turn on or off two log norm fit instead of one log norm fit
two_gaus = True
