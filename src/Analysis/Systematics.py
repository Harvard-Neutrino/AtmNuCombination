import numpy as np

####################
# Atmospheric flux #
####################

def FluxNormalization(x, experiment):
	return np.ones(experiment.NumberOfEvents) * x

def FluxNormalization_Below1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue<1] = nev[experiment.ETrue<1] * x
	return nev

def FluxNormalization_Above1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue>1] = nev[experiment.ETrue>1] * x
	return nev

def FluxTilt(x, experiment):
	E0Gam = 10 # GeV
	return (experiment.ETrue / E0Gam)**x

def NuNuBarRatio(x, experiment):
	nnbar = np.ones(experiment.NumberOfEvents)
	nnbar = nnbar[experiment.nuPDG<0] * x
	return nnbar

def FlavorRatio(x,experiment):
	eovermu = np.ones(experiment.NumberOfEvents)
	eovermu = eovermu[abs(experiment.nuPDG)==12] * x
	return eovermu

def ZenithFluxUp(x, experiment):
	zenith = np.ones(experiment.NumberOfEvents) 
	zenith = zenith[experiment.CosZTrue>=0] - x * np.tanh(experiment.CosZTrue[experiment.CosZTrue>=0])**2
	return zenith

def ZenithFluxDown(x, experiment):
	zenith = np.ones(experiment.NumberOfEvents) 
	zenith = zenith[experiment.CosZTrue<0] - x * np.tanh(experiment.CosZTrue[experiment.CosZTrue<0])**2
	return zenith

