import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQuIDS as nsq
import nuflux
import seaborn as sns
import scipy as scp

import util
from params import *
import classdef as cl

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})

def plot_contour(savename = "chi_sq_contour"):
	sin2t23 = np.sin(t23l) ** 2
	m31 = m31l

	X, Y = np.meshgrid(sin2t23, m31)
	Z = util.read_output()

	ogt23step = 0.001 * np.pi
	ogm31step = 0.01e-3
	Z2 = util.read_output(dir_name = "0117_IC", vt23l = np.arange(t23min, t23max + ogt23step, ogt23step), \
									vm31l = np.arange(m31min, m31max + ogm31step, ogm31step))

	Zinterp = scp.interpolate.interp2d(sin2t23, m31, Z, kind = 'cubic')
	newsin = np.arange(t23min, t23max + ogt23step, ogt23step)
	newx = np.sin(newsin) ** 2
	newy = np.arange(m31min, m31max + ogm31step, ogm31step)
	newz = Zinterp(newx, newy)
	XX, YY = np.meshgrid(newx, newy)

	fig, ax  = plt.subplots(figsize=(7,6))
	fig.suptitle("Chi-Sq Contour NH")
	ax.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax.set_ylabel(r"$m^2_{31}$")
	axim = ax.contour(XX,YY,newz,levels=[4.605, 5.991, 9.21],cmap=plt.cm.jet, linewidths = 0.8)
	axim = ax.contour(XX,YY,Z2,levels=[4.605, 5.991, 9.21], linewidths = 0.8)
	cb   = fig.colorbar(axim)
	fig.savefig("{}.png".format(savename), bbox_inches="tight")

plot_contour("0117_ORCA_coarse_interpolate_vs_IC.png")