import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQuIDS as nsq
import nuflux
import seaborn as sns

import util
from params import *

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})
def plot_contour(savename = "chi_sq_contour"):
	sin2t23 = np.sin(t23l) ** 2
	m31 = m31l

	X, Y = np.meshgrid(sin2t23, m31)
	Z = util.read_output()
	fig, ax  = plt.subplots(figsize=(7,6))
	fig.suptitle("Chi-Sq Contour NH")
	ax.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax.set_ylabel(r"$m^2_{31}$")
	axim = ax.contour(X,Y,Z,levels=[4.605, 5.991, 9.21],cmap=plt.cm.jet)
	cb   = fig.colorbar(axim)
	fig.savefig("{}.png".format(savename), bbox_inches="tight")

def plot_profile(var, savename = "profile"):
	fig2, ax2 = plt.subplots(figsize=(7,6))
	if var == 0: #theta
		profile = util.read_output()[0]
		x = np.sin(t23l) ** 2
		only_norm = util.read_output(dir_name = "1109_theta_profile_norm_min")[0]
		no_syst = util.read_output(dir_name = "../NoSystematics/1109_no_sys_the_profile")[0]
		fig2.suptitle(r"$\theta_{23} \chi^2$ profile")
		ax2.set_xlabel(r"$\sin^2{\theta_{23}}$")
		ax2.set_xlim([0.35, 0.65])
		ax2.set_ylim([0, 15])
	if var == 1: #m
		profile = util.read_output().T[0]
		only_norm = util.read_output(dir_name = "./1109_theta_profile_norm_min/").T[0]
		no_syst = util.read_output(dir_name = "../NoSystematics/1109_no_sys_m_profile").T[0]
		x = m31l
		fig2.suptitle(r"$\Delta m_{31}^2 \chi^2$ profile")
		ax2.set_xlabel(r"$\Delta m_{31}^2$")
		ax2.set_xlim([0.0024, 0.0026])
		ax2.set_ylim([0, 15])
	y = profile
	y_no_sys = no_syst
	y_only_norm = only_norm
	
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	ax2.plot(x, y, color ="green", label = "Norm + delta")
	ax2.plot(x, y_no_sys, '-.', color ="red", label = "No Syst")
	ax2.plot(x, y_only_norm, '--', color ="orange", label = "Only Norm")
	ax2.grid(True)
	ax2.legend()
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')

plot_profile(0, "1120_theta_norm+delta_min")