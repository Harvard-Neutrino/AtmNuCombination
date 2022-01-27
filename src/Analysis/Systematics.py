import numpy as np
import matplotlib.pyplot as plt

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
	return np.ones(experiment.NumberOfEvents) * (experiment.ETrue / E0Gam)**x

def NuNuBarRatio(x, experiment):
	nnbar = np.ones(experiment.NumberOfEvents)
	nnbar[experiment.nuPDG<0] = nnbar[experiment.nuPDG<0] * x
	return nnbar

def FlavorRatio(x,experiment):
	eovermu = np.ones(experiment.NumberOfEvents)
	eovermu[abs(experiment.nuPDG)==12] = eovermu[abs(experiment.nuPDG)==12] * x
	return eovermu

def ZenithFluxUp(x, experiment):
	zenith = np.ones(experiment.NumberOfEvents) 
	zenith[experiment.CosZTrue>=0] = zenith[experiment.CosZTrue>=0] - x * np.tanh(experiment.CosZTrue[experiment.CosZTrue>=0])**2
	return zenith

def ZenithFluxDown(x, experiment):
	zenith = np.ones(experiment.NumberOfEvents) 
	zenith[experiment.CosZTrue<0] = zenith[experiment.CosZTrue<0] - x * np.tanh(experiment.CosZTrue[experiment.CosZTrue<0])**2
	return zenith

def SKEnergyScale(x, experiment):
	modEReco = experiment.EReco * x
	w = experiment.Weight
	modW = np.ones(experiment.NumberOfEvents)
	for s in range(experiment.NumberOfSamples):
		bins = experiment.EnergyBins[s]
		for bn in range(bins.size -1):
			cond0 = (experiment.Sample==s) * (experiment.EReco>bins[bn]) * (experiment.EReco<bins[bn+1])
			cond1 = (experiment.Sample==s) * (modEReco>bins[bn]) * (modEReco<bins[bn+1])
			w1 = np.sum(w[cond1])
			w0 = np.sum(w[cond0])
			if w1==0:
				diff = 0
				modW[cond0] = modW[cond0]
			elif w0==0:
				diff = 1
			else:
				diff = w1 / w0
			modW[cond0] = modW[cond0] * diff
	return modW

