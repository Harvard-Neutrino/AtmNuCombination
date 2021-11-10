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
	ax.set_xlim([0.33, 0.67])
	ax.set_ylim([0.0022, 0.0028])
	axim = ax.contour(X,Y,Z,levels=[4.605, 5.991, 9.21],cmap=plt.cm.jet)
	cb   = fig.colorbar(axim)
	fig.savefig("{}.png".format(savename), bbox_inches="tight")

def plot_profile(var, savename = "profile"):
	profile = util.read_output()[0]
	if var == 0: #theta
		x = np.sin(t23l) ** 2
	if var == 1: #m
		x = m31l
	y = profile
	fig2, ax2 = plt.subplots(figsize=(7,6))
	fig2.suptitle(r"$\theta_{23} \chi^2$ profile (raw)")
	ax2.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax2.set_ylabel(r"$\chi^2_{NH}$")
	ax2.set_xlim([0.35, 0.65])
	ax2.set_ylim([0, 15])
	ax2.plot(x, y, color ="green")
	ax2.grid(True)
	fig2.savefig("{}.png".format(savename), bbox_inches='tight')
