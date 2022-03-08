import numpy as np

####################
#  Cross-section   #
####################

def XSecNuTau(x, experiment):
	tau = np.ones(experiment.NumberOfEvents)
	tau[np.abs(experiment.nuPDG)==16] *= x
	return tau 

def NCoverCC(x, experiment):
	nc = np.ones(experiment.NumberOfEvents)
	nc[experiment.CC==0] *= x 
	return nc

def NCHad(x, experiment):
	nc = np.ones(experiment.NumberOfEvents)
	nc[experiment.CC==0] *= x 
	return nc

def DIS(x,experiment):
	dis = np.ones(experiment.NumberOfEvents)
	cond = np.abs(experiment.Mode)>25 * experiment.CC
	dis[cond] *= x
	return dis

def CCQE(x,experiment):
	ccqe = np.ones(experiment.NumberOfEvents)
	ccqe[np.abs(experiment.Mode)==1] *= x
	return ccqe 

def CCQENuBarNu(x,experiment):
	ccqe = np.ones(experiment.NumberOfEvents)
	ccqe[experiment.Mode==-1] *= x
	return ccqe

def CCQEMuE(x,experiment):
	ccqe = np.ones(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)==1) * (np.abs(experiment.nuPDG)==14)
	ccqe[cond] *= x
	return ccqe

def CC1Pi_Pi0Pi(x,experiment):
	ccpi = np.ones(experiment.NumberOfEvents)
	ccpi[np.abs(experiment.Mode)==12] *= x
	return ccpi

def CC1Pi_NuBarNuE(x,experiment):
	ccpi = np.ones(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)>10) * (np.abs(experiment.Mode)<17) * (experiment.nuPDG==-12)
	ccpi[cond] *= x
	return ccpi

def CC1Pi_NuBarNuMu(x,experiment):
	ccpi = np.ones(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)>10) * (np.abs(experiment.Mode)<17) * (experiment.nuPDG==-14)
	ccpi[cond] *= x
	return ccpi

def CC1PiProduction(x,experiment):
	ccpi = np.ones(experiment.NumberOfEvents)
	cond = (np.abs(experiment.Mode)>10) * (np.abs(experiment.Mode)<17)
	ccpi[cond] *= x
	return ccpi

def CohPiProduction(x,experiment):
	ccpi = np.ones(experiment.NumberOfEvents)
	ccpi[np.abs(experiment.Mode)==16] *= x
	return ccpi
