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
res_e_true, res_e_reco, res_zen_true, res_zen_reco, res_W_mc, topology, pdg, \
									current_type, int_type = fakegen.generate()

if write:
	df_array = np.ndarray((17, len(res_e_true)))
	df_array[0] = [int(i) for i in pdg]
	df_array[1] = [int(i) for i in topology]
	df_array[2] = res_e_true
	df_array[3] = res_zen_true
	df_array[4] = input_MC["true_azimuth"]
	df_array[5] = res_e_reco
	df_array[6] = res_zen_reco
	df_array[7] = input_MC["reco_azimuth"]
	df_array[8] = [int(i) for i in int_type]
	df_array[9] = [int(i) for i in current_type]
	df_array[10] = res_W_mc
	df_array[11] = input_MC["xsec"]
	df_array[12] = input_MC["dxsec"]
	df_array[13] = input_MC["x"]
	df_array[14] = input_MC["y"]
	df_array[15] = input_MC["W"]
	df_array[16] = input_MC["Q2"]
	df = pd.DataFrame(df_array.T, columns = ["pdg", "pid", "true_energy", "true_zenith", \
		"true_azimuth", "reco_energy", "reco_zenith", "reco_azimuth", "interaction_type",\
		"current_type", "weight", "xsec", "dxsec", "x", "y", "W", "Q2"])
	df["pdg"].astype('int64')
	df["pid"].astype('int64')
	df["interaction_type"].astype('int64')
	df["current_type"].astype('int64')
	# df.to_csv("ORCA_pessimistic_Ereco.csv")
	df.to_csv("ORCA.csv")
