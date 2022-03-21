import numpy as np
import digitalizer as dgt
import MCGenerator as gen
import pandas as pd
from matplotlib import pyplot as plt
from params import *
import util

write = True

tracks = dgt.Digitalizer(input_track, input_scale)
cascades = dgt.Digitalizer(input_cascade, input_scale)

tracks.set_palette(0, -3, 100)
cascades.set_palette(0, -3, 100)

tracks.digitalize(22, 22)
cascades.digitalize(22, 22)

D_tracks = tracks.extracted
D_tracks[11][14] = 0.1 # only manual hardcode part
D_cascades = cascades.extracted

tracks.fit()
cascades.fit()

fakegen = gen.Generator(input_MC, tracks.gaussians, cascades.gaussians, x_bins)
res_e_true, res_e_reco, res_zen_true, res_zen_reco, res_W_mc, topology, pdg, int_type = fakegen.generate()

if write:
	df_array = np.ndarray((8, len(res_e_true)))
	df_array[0] = topology
	df_array[1] = pdg
	df_array[2] = res_e_true
	df_array[3] = res_e_reco
	df_array[4] = res_zen_true
	df_array[5] = res_zen_reco
	df_array[6] = res_W_mc
	df_array[7] = int_type
	df = pd.DataFrame(df_array.T, columns = ["pid", "pdg", "true_energy", "reco_energy", "true_zenith", "reco_zenith", "weight", "interaction_type"])
	df.to_csv("ORCA.csv")