import numpy as np 
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import scipy as scp
import scipy.interpolate

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})

# first extract the data
IO = pd.read_csv("./data/XiDCP_IO.dat", sep = " ")
NO = pd.read_csv("./data/XiDCP_NO.dat", sep = " ")

NO_nosys = pd.read_csv("./data/IC+SK+SKGd_DCP_Nosyst.dat", sep = " ")

# now plot them
def plot_dcp():
	DCP = IO["DCP"]
	XiNO = NO["Xi"]
	XiIO = IO["Xi"]
	XiNO_nosys = NO_nosys["X2"]
	DCP_nosys = NO_nosys["dCP"]

	fIO = scipy.interpolate.interp1d(DCP, XiIO)
	fNO = scipy.interpolate.interp1d(DCP, XiNO)
	fNO_nosys = scipy.interpolate.interp1d(DCP_nosys, XiNO_nosys)

	DCPnew = np.arange(0, 2 * np.pi, 0.01)
	IOnew = fIO(DCPnew)
	NOnew = fNO(DCPnew)
	NO_nosys_new = fNO_nosys(DCPnew)

	fig, ax  = plt.subplots(figsize=(10,8))
	fig.suptitle(r"$\delta_{CP}$ Sensitivity: SK+SKGd(5yrs)+ICUp(5yrs) (Assuming NO, preliminary)")
	ax.set_xlabel(r"$\delta_{CP}$")
	ax.set_ylabel(r"$\Delta X^2$")
	ax.plot(DCPnew, IOnew, label = "IO", color = "green")
	ax.plot(DCPnew, NOnew, label = "NO", color = "red")
	ax.plot(DCPnew, NO_nosys_new, label = "NO, no syst", color = "red", linestyle = "--")
	ax.axhline(y=1, color='black', linestyle='--', linewidth = 1.0)
	ax.axhline(y=4, color='black', linestyle='--', linewidth = 1.0)
	ax.axhline(y=9, color='black', linestyle='--', linewidth = 1.0)
	ax.axhline(y=16, color='black', linestyle='--', linewidth = 1.0)
	ax.axhline(y=25, color='black', linestyle='--', linewidth = 1.0)
	ax.legend()
	# plt.show()
	fig.savefig("DCP", bbox_inches="tight")
	plt.close()

plot_dcp()