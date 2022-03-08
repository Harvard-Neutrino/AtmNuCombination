import numpy as np

####################
# Atmospheric flux #
####################

def Diff_FluxNormalization(x, experiment):
	return 1

def Diff_FluxNormalization_Below1GeV(x, experiment):
	nev = np.zeros(experiment.NumberOfEvents)
	nev[experiment.ETrue<1] = 1
	return nev

def Diff_FluxNormalization_Above1GeV(x, experiment):
	nev = np.zeros(experiment.NumberOfEvents)
	nev[experiment.ETrue>1] = 1
	return nev

def Diff_FluxTilt(x, experiment):
	E0Gam = 10 # GeV
	return np.ones(experiment.NumberOfEvents) * (experiment.ETrue / E0Gam)**x * np.log(experiment.ETrue / E0Gam)

def Diff_NuNuBarRatio(x, experiment):
	nnbar = np.zeros(experiment.NumberOfEvents)
	nnbar[experiment.nuPDG<0] = 1
	return nnbar

def Diff_FlavorRatio(x,experiment):
	eovermu = np.zeros(experiment.NumberOfEvents)
	eovermu[abs(experiment.nuPDG)==12] = 1
	return eovermu

def Diff_ZenithFluxUp(x, experiment):
	zenith = np.zeros(experiment.NumberOfEvents) 
	zenith[experiment.CosZTrue>=0] = - np.tanh(experiment.CosZTrue[experiment.CosZTrue>=0])**2
	return zenith

def Diff_ZenithFluxDown(x, experiment):
	zenith = np.zeros(experiment.NumberOfEvents) 
	zenith[experiment.CosZTrue<0] = - np.tanh(experiment.CosZTrue[experiment.CosZTrue<0])**2
	return zenith