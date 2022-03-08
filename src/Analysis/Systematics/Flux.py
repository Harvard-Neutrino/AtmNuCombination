import numpy as np

####################
# Atmospheric flux #
####################

def FluxNormalization(x, experiment):
	return x

def FluxNormalization_Below1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue<1] = x
	return nev

def FluxNormalization_Above1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue>1] = x
	return nev

def FluxTilt(x, experiment):
	E0Gam = 10 # GeV
	return np.ones(experiment.NumberOfEvents) * (experiment.ETrue / E0Gam)**x

def NuNuBarRatio(x, experiment):
	nnbar = np.ones(experiment.NumberOfEvents)
	nnbar[experiment.nuPDG<0] = x
	return nnbar

def FlavorRatio(x,experiment):
	eovermu = np.ones(experiment.NumberOfEvents)
	eovermu[abs(experiment.nuPDG)==12] = x
	return eovermu

def ZenithFluxUp(x, experiment):
	zenith = np.ones(experiment.NumberOfEvents) 
	zenith[experiment.CosZTrue>=0] = zenith[experiment.CosZTrue>=0] - x * np.tanh(experiment.CosZTrue[experiment.CosZTrue>=0])**2
	return zenith

def ZenithFluxDown(x, experiment):
	zenith = np.ones(experiment.NumberOfEvents) 
	zenith[experiment.CosZTrue<0] = zenith[experiment.CosZTrue<0] - x * np.tanh(experiment.CosZTrue[experiment.CosZTrue<0])**2
	return zenith