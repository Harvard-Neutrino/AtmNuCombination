import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit
import math
# from Effective import ICEffectiveAnalysis as ICEff
# from Effective import ORCAEffectiveAnalysis as ORCAEff
import digitalizer as dgt
import util
from util import gaussian
# import analyze as anl
from params import *

input_file = pd.read_csv(savename)
original = pd.read_csv("neutrino_mc.csv")
pid = input_file["pid"]
e_true = input_file["true_energy"]
e_reco = input_file["reco_energy"]
zen_true = input_file["true_zenith"]
zen_reco = input_file["reco_zenith"]
weights = input_file["weight"]

def check_int_weight():
	cnt0 = 0
	cnt1 = 0
	cnt2 = 0
	for i in range(len(weights)):
	    if pid[i] == 2 and weights[i] != 0:
	        print(i)
	    if pid[i] == 0 and weights[i] == 0:
	        print(i)
	    if pid[i] == 0:
	        cnt0 += 1
	    if pid[i] == 1:
	        cnt1 += 1
	    if pid[i] == 2:
	        cnt2 += 1
	print(cnt0, cnt1, cnt2)

def select_ICMC():
	mask = []
	for i in range(len(original["weight"])):
		if original["true_energy"][i] >= 54 or original["true_energy"][i] <= 1:
			mask.append(False)
		else:
			mask.append(True)

	input_MC = original
	df_array = np.ndarray((17, len(input_MC["pdg"][mask])))
	df_array[0] = input_MC["pdg"][mask]
	df_array[1] = input_MC["pid"][mask]
	df_array[2] = input_MC["true_energy"][mask]
	df_array[3] = input_MC["true_zenith"][mask]
	df_array[4] = input_MC["true_azimuth"][mask]
	df_array[5] = input_MC["reco_energy"][mask]
	df_array[6] = input_MC["reco_zenith"][mask]
	df_array[7] = input_MC["reco_azimuth"][mask]
	df_array[8] = input_MC["interaction_type"][mask]
	df_array[9] = input_MC["current_type"][mask]
	df_array[10] = input_MC["weight"][mask]
	df_array[11] = input_MC["xsec"][mask]
	df_array[12] = input_MC["dxsec"][mask]
	df_array[13] = input_MC["x"][mask]
	df_array[14] = input_MC["y"][mask]
	df_array[15] = input_MC["W"][mask]
	df_array[16] = input_MC["Q2"][mask]
	df = pd.DataFrame(df_array.T, columns = ["pdg", "pid", "true_energy", "true_zenith", \
		"true_azimuth", "reco_energy", "reco_zenith", "reco_azimuth", "interaction_type",\
		"current_type", "weight", "xsec", "dxsec", "x", "y", "W", "Q2"])


	df.to_csv("selected_ICMC.csv")

def remove_morph(filename, n = 2):
	infile = pd.read_csv(filename)
	new_df = infile.copy()
	for i in range(len(new_df["pid"])):
		if new_df["pid"][i] == n:
			new_df["weight"][i] = 0
	new_df.to_csv("{}_without_int.csv".format(filename))

# remove_morph(savename)

def apply_weight(source, target, save):
	from_file = pd.read_csv(source)
	to_file = pd.read_csv(target)
	new_df = to_file.copy()
	for i in range(len(to_file["weight"])):
		if new_df["weight"][i] != 0:
			new_df["weight"][i] = from_file["weight"][i]
	new_df.to_csv("{}.csv".format(save))

# apply_weight("0623ORCA_morph_test.csv","0624ORCA_morph_wo_rw_without_int.csv", "0627ORCA_applied_rw_to_morph_wo_int")
# apply_weight("0627ORCA_no_morph_lin_fit.csv", "0623ORCA_no_morph_control.csv", "0627ORCA_no_morph_control_applied_lin_fit_rw")

def do_nothing():
	return

# 0628 checks what's the difference between two supposedly similar files
def check_why():
	no_morph = pd.read_csv("0623ORCA_no_morph_control.csv")
	yes_morph = pd.read_csv("0628ORCA_applied_morph_to_control_without_int.csv")
	# first check total number of entries:
	if len(no_morph["weight"]) == len(yes_morph["weight"]):
		# print("check 1: yay!")
		do_nothing()
	else:
		print("check 1: oh no!")

	# then check tracks and cascades stay the same
	for i in range(len(no_morph["pid"])):
		if yes_morph["pid"][i] == 2:
			if yes_morph["weight"][i] == 0:
				# print("check 2.1: yay!")
				do_nothing()
			else:
				print("check 2.1: no!")
		else:
			if yes_morph["pid"][i] == no_morph["pid"][i]:
				# print("check 2.2: yay!")
				do_nothing()
			else:
				print("check 2.2: no!")

	# then check all the weight are the same:
	for i in range(len(no_morph["weight"])):
		if no_morph["weight"][i] - yes_morph["weight"][i] <= 0.0001 * no_morph["weight"][i]:
			# print("check 3: yay!")
			do_nothing()
		else:
			if yes_morph["weight"][i] != 0:
				print(i)
				print(no_morph["weight"][i])
				print(yes_morph["weight"][i])

	# then check all the reco energies are the same
	for i in range(len(no_morph["true_energy"])):
		if no_morph["true_energy"][i] - yes_morph["true_energy"][i] <= 0.0001 * no_morph["true_energy"][i]:
			# print("check 3: yay!")
			do_nothing()
		else:
			print(i)

# looks like the two files passed the check ???
check_why()

def apply_morph(source, target, save):
	from_file = pd.read_csv(source)
	to_file = pd.read_csv(target)
	new_df = to_file.copy()
	for i in range(len(to_file["pid"])):
		new_df["pid"][i] = from_file["pid"][i]
	new_df.to_csv("{}.csv".format(save))

# apply_morph("0627ORCA_applied_rw_to_morph_wo_int.csv", "0623ORCA_no_morph_control.csv", "0628ORCA_applied_morph_to_control")
# remove_morph("0628ORCA_applied_morph_to_control.csv")

