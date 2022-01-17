import numpy as np 
import pandas as pd

input_scale = "./ORCA/scale.png"
input_track = "./ORCA/migrationTracks.png"
input_cascade = "./ORCA/migrationCascades.png"

input_MC = pd.read_csv("neutrino_mc.csv")

x_bins = np.logspace(np.log10(1.85), np.log10(53), 23)