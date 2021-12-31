import numpy as np

def FluxNormalization(x, experiment):
	return np.ones(experiment.NumberOfEvents) * x

def FluxNormalization_Below1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue<1] = nev[experiment.ETrue<1] * x
	return nev

def SKFluxNormalization_Above1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue>1] = nev[experiment.ETrue>1] * x
	return nev

def SKFluxNormalization_Below1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue<1] = nev[experiment.ETrue<1] * x
	return nev

def ICFluxNormalization_Above1GeV(x, experiment):
	nev = np.ones(experiment.NumberOfEvents)
	nev[experiment.ETrue>1] = nev[experiment.ETrue>1] * x
	return nev