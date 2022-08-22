import numpy as np
import pandas as pd
import util
import sys
import math
from params import *

# this file multiplies each row of a given csv NUM times
NUM = 2
OUTPUT_NAME = "{}xneutrino_mc.csv".format(NUM)

input_mc = pd.read_csv("neutrino_mc.csv")
length = len(input_mc["pdg"])
tot_count = length
new_array = np.ndarray((17, NUM * length))

for i in range(length):
	cur_count = 0 
	new_array[0][i] = input_mc["pdg"][i]
	new_array[1][i] = input_mc["pid"][i]
	new_array[2][i] = input_mc["true_energy"][i]
	new_array[3][i] = input_mc["true_zenith"][i]
	new_array[4][i] = input_mc["true_azimuth"][i]
	new_array[5][i] = input_mc["reco_energy"][i]
	new_array[6][i] = input_mc["reco_zenith"][i]
	new_array[7][i] = input_mc["reco_azimuth"][i]
	new_array[8][i] = input_mc["interaction_type"][i]
	new_array[9][i] = input_mc["current_type"][i]
	new_array[10][i] = input_mc["weight"][i] / NUM
	new_array[11][i] = input_mc["xsec"][i]
	new_array[12][i] = input_mc["dxsec"][i]
	new_array[13][i] = input_mc["x"][i]
	new_array[14][i] = input_mc["y"][i]
	new_array[15][i] = input_mc["W"][i]
	new_array[16][i] = input_mc["Q2"][i]
	for j in range(NUM - 1):
		new_array[0][tot_count + j] = input_mc["pdg"][i]
		new_array[1][tot_count + j] = input_mc["pid"][i]
		new_array[2][tot_count + j] = input_mc["true_energy"][i]
		new_array[3][tot_count + j] = input_mc["true_zenith"][i]
		new_array[4][tot_count + j] = input_mc["true_azimuth"][i]
		new_array[5][tot_count + j] = input_mc["reco_energy"][i]
		new_array[6][tot_count + j] = input_mc["reco_zenith"][i]
		new_array[7][tot_count + j] = input_mc["reco_azimuth"][i]
		new_array[8][tot_count + j] = input_mc["interaction_type"][i]
		new_array[9][tot_count + j] = input_mc["current_type"][i]
		new_array[10][tot_count + j] = input_mc["weight"][i] / NUM
		new_array[11][tot_count + j] = input_mc["xsec"][i]
		new_array[12][tot_count + j] = input_mc["dxsec"][i]
		new_array[13][tot_count + j] = input_mc["x"][i]
		new_array[14][tot_count + j] = input_mc["y"][i]
		new_array[15][tot_count + j] = input_mc["W"][i]
		new_array[16][tot_count + j] = input_mc["Q2"][i]
	tot_count += NUM - 1

print("the original mc contained {} total number of events".format(length))
print("the new multiplied mc now contains {} total number of events".format(tot_count))

df = pd.DataFrame(new_array.T, columns = ["pdg", "pid", "true_energy", "true_zenith", \
		"true_azimuth", "reco_energy", "reco_zenith", "reco_azimuth", "interaction_type",\
		"current_type", "weight", "xsec", "dxsec", "x", "y", "W", "Q2"])
df.to_csv(OUTPUT_NAME)

