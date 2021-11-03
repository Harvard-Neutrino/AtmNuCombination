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
