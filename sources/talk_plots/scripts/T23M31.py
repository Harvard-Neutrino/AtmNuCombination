import numpy as np 
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import scipy as scp
import scipy.interpolate

matplotlib.rcParams.update({'font.size': 23})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})

# first extract the data
with_sys = pd.read_csv("./../data/XiDM31_T23_Sys.dat", sep = " ")
no_sys = pd.read_csv("./../data/XiDM31_T23_Nosys.dat", sep = " ")

def plot_t23_dm31(incsys):
	if incsys:
		sys = with_sys
		plot_name = "Chi-Sq Contour (NO, with systematic))"
		savename = "./../plots/T23DM31_sys"
	else:
		sys = no_sys
		plot_name = "Chi-Sq Contour (NO, without systematic))"
		savename = "./../plots/T23DM31_nosys"
	theta23 = sys["Theta23"]
	dm31 = sys["DM31"]
	chisq = sys["ChiSq"]
	
	X, Y = np.meshgrid(theta23, dm31)
	Z = chisq

	t23min = 0.3
	t23max = 0.7
	m31min = 0.00235
	m31max = 0.00265
	ogt23step = 0.001 * np.pi
	ogm31step = 0.01e-3

	Zinterp = scp.interpolate.interp2d(theta23, dm31, Z, kind = 'cubic')
	newx = np.arange(t23min, t23max + ogt23step, ogt23step)
	newy = np.arange(m31min, m31max + ogm31step, ogm31step)
	newz = Zinterp(newx, newy)
	XX, YY = np.meshgrid(newx, newy)

	fig, ax  = plt.subplots(figsize=(13,11))
	fig.suptitle(plot_name)
	ax.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax.set_ylabel(r"$\Delta m^2_{31}$")
	axim = ax.contour(XX,YY,newz,levels=[4.605, 5.991, 9.21],cmap=plt.cm.jet, linewidths = 1.2)
	cb   = fig.colorbar(axim)
	# plt.show()
	fig.savefig(savename, bbox_inches="tight")
	plt.close()

def plot_t23_dm31_together():

	plot_name = "Sensitivity w/ and w/o Systematics"
	savename = "./../plots/T23DM31_together"

	sys1 = with_sys
	sys2 = no_sys

	theta23 = sys1["Theta23"]
	dm31 = sys1["DM31"]
	chisq1 = sys1["ChiSq"]

	theta232 = sys2["Theta23"]
	dm312 = sys2["DM31"]
	chisq2 = sys2["ChiSq"]
	
	X, Y = np.meshgrid(theta23, dm31)
	Z1 = chisq1
	Z2 = chisq2

	t23min = 0.3
	t23max = 0.7
	m31min = 0.00235
	m31max = 0.00265
	ogt23step = 0.001 * np.pi
	ogm31step = 0.01e-3

	Z1interp = scp.interpolate.interp2d(theta23, dm31, Z1, kind = 'cubic')
	Z2interp = scp.interpolate.interp2d(theta232, dm312, Z2, kind = 'cubic')
	newx = np.arange(t23min, t23max + ogt23step, ogt23step)
	newy = np.arange(m31min, m31max + ogm31step, ogm31step)
	newz1 = Z1interp(newx, newy)
	newz2 = Z2interp(newx, newy)
	XX, YY = np.meshgrid(newx, newy)

	fig, ax  = plt.subplots(figsize=(11,8))
	fig.suptitle(plot_name)
	ax.set_xlabel(r"$\sin^2{\theta_{23}}$")
	ax.set_ylabel(r"$\Delta m^2_{31}$ [eV^2]")
	colors = ['sienna', 'lawngreen', 'olivedrab']
	axim = ax.contour(XX,YY,newz1,levels=[4.605, 5.991, 9.21],colors = colors, \
									linewidths = 1.2, label = "with systematics")
	axim2 = ax.contour(XX,YY,newz2, linestyles = "dashed", levels=[4.605, 5.991, 9.21], \
										colors = colors, linewidths = 1.2, label = "without systematics")
	h1,_ = axim.legend_elements()
	ax.legend([h1[0], h1[1], h1[2]], [r'SK+SKGd(5yrs)+ICUp(5yrs) @ 90% CL',\
	r'SK+SKGd(5yrs)+ICUp(5yrs) @ 95% CL',  r'SK+SKGd(5yrs)+ICUp(5yrs) @ 99% CL'])
	ax.ticklabel_format(style = 'sci', scilimits = (-2,1))
	# ax.legend()
	# plt.show()
	fig.savefig(savename, bbox_inches="tight")
	plt.close()

# plot_t23_dm31(True)
# plot_t23_dm31(False)
plot_t23_dm31_together()
