import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQuIDS as nsq
import nuflux
import seaborn as sns

import sensitivity as sst
from params import *

def hist_truth():
	rate_weight_truth, energy_hist_truth, energy_bins_truth = sst.get_rated_weight_truth(gamma = -1)
	input_data["rate_weight"] = rate_weight_truth
	
	# Plot the energy distribution
	print("plotting energy distributions")
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("Reco Energy Rated Distribution")
	ax.hist(input_data["reco_energy"][cascade_mask], bins=E_bin_plot, \
		  weights=input_data["rate_weight"][cascade_mask], \
		  label=r"$\nu_{All, Cascade}$", color="blue", histtype="step")
	ax.hist(input_data["reco_energy"][track_mask], bins=E_bin_plot, \
		  weights=input_data["rate_weight"][track_mask], \
		  label=r"$\nu_{All, Track}$", color="red", histtype="step")
	ax.set_xscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{reco}}$ [GeV]")
	ax.set_xlim(1, 1000)
	ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
					 useOffset=None, useLocale=None, useMathText=None)
	ax.set_ylabel("Rate [3 Years]")
	ax.grid(True)
	ax.legend()
	# fig.show()
	fig.savefig("TruthFit(norm = 1, gamma = -1).png", bbox_inches='tight')
	
	#fig.savefig("Hist(Truth, no nuisance).png", bbox_inches='tight')

def hist_nuisance(dm, th, norm, nudelta):
	rate_weight = sst.get_energy_bins(dm, th, top = 2, norm = norm, delta = nudelta)
	input_data["rate_weight"] = rate_weight
	
	# Plot the energy distribution
	print("plotting energy distributions")
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("DM = {}, s^2th = {}, norm = {}, delta = {}".format(dm, th, norm, nudelta))
	ax.hist(input_data["reco_energy"][cascade_mask], bins=E_bin_plot, \
		  weights=input_data["rate_weight"][cascade_mask], \
		  label=r"$\nu_{All, Cascade}$", color="blue", histtype="step")
	ax.hist(input_data["reco_energy"][track_mask], bins=E_bin_plot, \
		  weights=input_data["rate_weight"][track_mask], \
		  label=r"$\nu_{All, Track}$", color="red", histtype="step")
	ax.set_xscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{reco}}$ [GeV]")
	ax.set_xlim(1, 1000)
	ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
					 useOffset=None, useLocale=None, useMathText=None)
	ax.set_ylabel("Rate [3 Years]")
	ax.grid(True)
	ax.legend()
	# fig.show()
	fig.savefig("Hist(DM = {}, s^2th = {}, norm = {}, delta = {}).png".format(dm, th, norm, nudelta), bbox_inches='tight')

def distribution_plots():
	print("plotting sanity check distributions")

	rate_weight_truth, energy_hist_truth, energy_bins_truth = sst.get_rated_weight_truth()
	input_data["rate_weight"] = rate_weight_truth
	
	# Plot the energy distribution
	print("plotting energy distributions")
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("Reco Energy Rated Distribution")
	ax.hist(input_data["reco_energy"][cascade_mask], bins=E_bin_plot, \
		  weights=input_data["rate_weight"][cascade_mask], \
		  label=r"$\nu_{All, Cascade}$", color="blue", histtype="step")
	ax.hist(input_data["reco_energy"][track_mask], bins=E_bin_plot, \
		  weights=input_data["rate_weight"][track_mask], \
		  label=r"$\nu_{All, Track}$", color="red", histtype="step")
	ax.set_xscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{reco}}$ [GeV]")
	ax.set_xlim(1, 1000)
	ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
					 useOffset=None, useLocale=None, useMathText=None)
	ax.set_ylabel("Rate [3 Years]")
	ax.grid(True)
	ax.legend()
	fig.savefig("Rate_For_Sensitivity.png", bbox_inches='tight')


	# Plot the angle distribution
	print("plotting zenith distributions")
	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("Reco Zenith Rated Distribution")
	ax.hist(np.cos(input_data["reco_zenith"][cascade_mask]), bins=cos_bin_plot, \
		  weights=input_data["rate_weight"][cascade_mask], \
		  label=r"$\nu_{All, Cascade}$", color="blue", histtype="step")
	ax.hist(np.cos(input_data["reco_zenith"][track_mask]), bins=cos_bin_plot, \
		  weights=input_data["rate_weight"][track_mask], \
		  label=r"$\nu_{All, Track}$", color="red", histtype="step")
	ax.set_xlabel(r"$\cos{\theta, \rm{reco}}$")
	ax.set_xlim(-1, 1)
	ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
					 useOffset=None, useLocale=None, useMathText=None)
	ax.set_ylabel("Rate [3 Years]")
	ax.grid(True)
	ax.legend()
	fig.savefig("Zenith_Rate_For_Sensitivity.png", bbox_inches='tight')

	# Plot the overall distribution
	print("plotting 2D distributions")
	counts, _, _ = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=[energy_bins_fine, cos_bin_plot], \
		  weights=input_data["rate_weight"])

	fig, ax = plt.subplots(figsize=(7,6))
	fig.suptitle("Reco Energy and Zenith Rated Distribution")
	ax.pcolormesh(energy_bins_fine, cos_bin_plot, counts.T)
	ax.set_xscale('log')
	ax.set_xlabel(r"$E_{\nu,\rm{reco}}$ [GeV]")
	ax.set_ylabel(r"$\cos{\theta, \rm{reco}}$")
	ax.set_xlim(1, 100)
	ax.set_ylim(-1, 1)
	ax.legend()
	fig.savefig("2D_Rate_For_Sensitivity.png", bbox_inches='tight')



# Plots contour of t23 and m31 chi-sq for NH
def plot_contour_chi(savename = "chi_sq_contour"):
	print("plotting chi squared contour")
	rate_weight_truth_0, energy_hist_truth_0, energy_bins_truth_0 = sst.get_rated_weight_truth(top = 0)
	rate_weight_truth_1, energy_hist_truth_1, energy_bins_truth_1 = sst.get_rated_weight_truth(top = 1)

	
	# # try some bigger values for the step
	# bigstep = True
	# if bigstep:
	# 	t23step = 0.02 * np.pi
	# 	m31step = 0.20e-3

	sin2t23 = np.arange(t23min, t23max + t23step, t23step)
	for i in range(len(sin2t23)):
		sin2t23[i] = np.sin(sin2t23[i]) ** 2
	m31 = np.arange(m31min, m31max + m31step, m31step)

	X, Y = np.meshgrid(sin2t23, m31)
	# Z = sst.get_chisq(X, Y, energy_hist_truth, top = top)

	Z2 = np.ndarray((len(m31), len(sin2t23)))
	for i in range(len(m31)):
		for j in range(len(sin2t23)):
			Z2[i][j] = sst.get_chisq(np.arcsin(np.sqrt(sin2t23[j])), m31[i],  energy_hist_truth_0, top = 0) \
						+ sst.get_chisq(np.arcsin(np.sqrt(sin2t23[j])), m31[i],  energy_hist_truth_1, top = 1)
			# print("m31, ", m31[i])
			# print("t23, ", np.arcsin(np.sqrt(sin2t23[j])))
			# print("sin2t23, ", sin2t23[j])
			# print("chisq, ", Z2[i][j])

	fig4, ax4  = plt.subplots(figsize=(7,6))
	fig4.suptitle("Chi-Sq Contour NH")
	ax4.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax4.set_ylabel(r"$m^2_{31}$")
	axim = ax4.contour(X,Y,Z2,levels=[4.605, 5.991, 9.21],cmap=plt.cm.jet)
	cb   = fig4.colorbar(axim)

	fig4.savefig("{}.png".format(savename), bbox_inches="tight")

# plot un-normalized chisq for NH, probing values of t23, not minimizing over m31
def plot_t23_chi_raw_profile(savename = "t23_chi_sq_profile_raw_new", top = 0):
	print("plotting t23 chi profile")
	x = np.sin(t23l) ** 2
	y = sst.get_t23_chi_profile(top = top)
	fig2, ax2 = plt.subplots(figsize=(7,6))
	fig2.suptitle(r"$\theta_{23} \chi^2$ profile (raw)")
	ax2.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	ax2.set_yscale("log")
	ax2.plot(x, y, color ="green")
	ax2.grid(True)
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')

# plot un-normalized chisq for NH, probing values of t23, minimizing over m31
def plot_t23_min_chi_profile(savename = "t23_chi_sq_min_profile", top = 0):
	print("plotting t23 min chi profile")
	x = np.sin(t23l) ** 2
	y = sst.get_t23_min_chi_profile(top = top)
	fig2, ax2 = plt.subplots(figsize=(7,6))
	fig2.suptitle(r"$\theta_{23} \chi^2$ profile (minimized)")
	ax2.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	ax2.set_yscale("log")
	ax2.plot(x, y, color ="green")
	ax2.grid(True)
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')

# plot un-normalized chisq for NH, probing values of t23, not minimizing over m31, but for all topologies
def plot_t23_chi_raw_profile_all_top(savename = "t23_chi_sq_profile_raw_all_top"):
	print("plotting t23 chi profile all topologies")
	x = np.sin(t23l) ** 2
	y0 = sst.get_t23_chi_profile(top = 0)
	y1 = sst.get_t23_chi_profile(top = 1)
	y2 = list()
	y2[:] = y0[:]
	for i in range(len(y2)):
		y2[i] = y0[i] + y1[i]
	fig2, ax2 = plt.subplots(figsize=(7,6))
	fig2.suptitle(r"$\theta_{23} \chi^2$ profile (raw)")
	ax2.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	# ax2.set_yscale("log")
	ax2.plot(x, y0, color ="green", label = "cascades")
	ax2.plot(x, y1, color ="red", label = "tracks")
	ax2.plot(x, y2, color ="blue", label = "all")
	ax2.grid(True)
	ax2.legend()
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')

# plot un-normalized chisq for NH, probing values of t23, not minimizing over m31, but for all topologies
def plot_t23_min_chi_profile_all_top(savename = "t23_min_chi_sq_profile_all_top"):
	print("plotting t23 chi profile all topologies")
	x = np.sin(t23l) ** 2
	y0 = sst.get_t23_min_chi_profile(top = 0)
	y1 = sst.get_t23_min_chi_profile(top = 1)
	y2 = sst.get_t23_min_chi_profile(top = 2)
	fig2, ax2 = plt.subplots(figsize=(7,6))
	fig2.suptitle(r"$\theta_{23} \chi^2$ profile (minimized)")
	ax2.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	ax2.set_yscale("log")
	ax2.plot(x, y0, color ="green", label = "cascades")
	ax2.plot(x, y1, color ="red", label = "tracks")
	ax2.plot(x, y2, color ="blue", label = "all")
	ax2.grid(True)
	ax2.legend()
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')

# plot un-normalized chisq for NH, probing values of t23, not minimizing over m31
def plot_m31_chi_raw_profile(savename = "m31_chi_sq_profile_raw_new", top = 0):
	print("plotting m31 chi profile")
	x = m31l
	y = sst.get_m31_chi_profile(top = top)
	fig2, ax2 = plt.subplots(figsize=(7,6))
	fig2.suptitle(r"$m_{31} \chi^2$ profile (raw)")
	ax2.set_xlabel(r"$m_{31}$")
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	ax2.set_yscale("log")
	ax2.plot(x, y, color ="green")
	ax2.grid(True)
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')

# plot un-normalized chisq for NH, probing values of m31, not minimizing over t23, but for all topologies
def plot_m31_chi_raw_profile_all_top(savename = "m31_chi_sq_profile_raw_all_top_new"):
	print("plotting m31 chi profile all topologies")
	x = m31l
	y0 = sst.get_m31_chi_profile(top = 0)
	y1 = sst.get_m31_chi_profile(top = 1)
	y2 = list()
	y2[:] = y0[:]
	for i in range(len(y2)):
		y2[i] = y0[i] + y1[i]
	fig2, ax2 = plt.subplots(figsize=(7,6))
	fig2.suptitle(r"$m_{31} \chi^2$ profile (raw)")
	ax2.set_xlabel(r"$m_{31}$")
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	ax2.set_yscale("log")
	ax2.plot(x, y0, color ="green", label = "cascades")
	ax2.plot(x, y1, color ="red", label = "tracks")
	ax2.plot(x, y2, color ="blue", label = "all")
	ax2.grid(True)
	ax2.legend()
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')

# Plot the resolution
def plot_resolution():
	# new_weight = np.zeros_like(input_data["weight"])
	# for i in range(len(input_data["weight"])):
	# 	if input_data["true_energy"][i] >= 10^(1.0) & input_data["true_energy"][i] < 10^(1.1):
	# 		new_weight[i] = 1
	# input_data["new_weight"] = new_weight
	fig, ax = plt.subplots(figsize = (7,6))
	fig.suptitle("Energy Resolution")
	ax.set_xlabel("Reco Energy")
	ax.set_ylabel("Normalized Rate")
	ax.hist(input_data["reco_energy"][reso_mask0], bins=energy_bins_fine, \
		  label=r"$\nu_{All, E_{True} \in (10^1, 10^{1.1}))]}$", color="blue", histtype="step")
	# sns.histplot(input_data, x = "reco_energy", element = "step", weights = "new_weight", ax = ax)
	fig.savefig("resolution.png", bbox_inches = "tight")
