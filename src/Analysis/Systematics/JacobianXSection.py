import numpy as np

####################
#  Cross-section   #
####################

def Diff_XSecNuTau(x, experiment):
	tau = np.zeros(experiment.NumberOfEvents)
	tau[np.abs(experiment.nuPDG)==16] = 1
	return tau 

def Diff_NCoverCC(x, experiment):
	nc = np.zeros(experiment.NumberOfEvents)
	nc[experiment.CC==0] = 1
	return nc

def Diff_NCHad(x, experiment):
	nc = np.zeros(experiment.NumberOfEvents)
	nc[experiment.CC==0] = 1 
	return nc

def Diff_DIS(x,experiment):
	dis = np.zeros(experiment.NumberOfEvents)
	cond = np.abs(experiment.Mode)>25 * experiment.CC
	dis[cond] = 1
	return dis

def Diff_CCQE(x,experiment):
	ccqe = np.zeros(experiment.NumberOfEvents)
	ccqe[np.abs(experiment.Mode)==1] = 1
	return ccqe 

def Diff_CCQENuBarNu(x,experiment):
	ccqe = np.zeros(experiment.NumberOfEvents)
	ccqe[experiment.Mode==-1] = 1
	return ccqe

def Diff_CCQEMuE(x,experiment):
	ccqe = np.zeros(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)==1) * (np.abs(experiment.nuPDG)==14)
	ccqe[cond] = 1
	return ccqe

def Diff_CC1Pi_Pi0Pi(x,experiment):
	ccpi = np.zeros(experiment.NumberOfEvents)
	ccpi[np.abs(experiment.Mode)==12] = 1
	return ccpi

def Diff_CC1Pi_NuBarNuE(x,experiment):
	ccpi = np.zeros(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)>10) * (np.abs(experiment.Mode)<17) * (experiment.nuPDG==-12)
	ccpi[cond] = 1
	return ccpi

def Diff_CC1Pi_NuBarNuMu(x,experiment):
	ccpi = np.zeros(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)>10) * (np.abs(experiment.Mode)<17) * (experiment.nuPDG==-14)
	ccpi[cond] = 1
	return ccpi

def Diff_CC1PiProduction(x,experiment):
	ccpi = np.zeros(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)>10) * (np.abs(experiment.Mode)<17)
	ccpi[cond] = 1
	return ccpi

def Diff_CohPiProduction(x,experiment):
	ccpi = np.zeros(experiment.NumberOfEvents)
	ccpi[np.abs(experiment.Mode)==16] = 1
	return ccpi

