import numpy as np 
import pandas as pd

input_scale = "./ORCA_Results/scale.png"
input_track = "./ORCA_Results/migrationTracks.png"
input_cascade = "./ORCA_Results/migrationCascades.png"

input_MC = pd.read_csv("neutrino_mc.csv")

x_bins = np.logspace(np.log10(1.85), np.log10(54), 23)