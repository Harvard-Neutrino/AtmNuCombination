import numpy as np

####################
####### ORCA #######
####################

def ORCAEnergyScale(x, experiment):
	# if np.abs(x-1)<5e-4: 
	# 	return np.zeros(experiment.NumberOfBins)
	return experiment.Exp_wBinIt(1,shift_E=x) / experiment.weightOscBF_binned - 1

def ShowerClassification(x,experiment):
	ev = np.ones(experiment.NumberOfEvents)
	ev[experiment.Sample==0] = x
	return experiment.Exp_wBinIt(ev) / experiment.weightOscBF_binned - 1

def IntermediateClassification(x,experiment):
	ev = np.ones(experiment.NumberOfEvents)
	ev[experiment.Sample==1] = x
	return experiment.Exp_wBinIt(ev) / experiment.weightOscBF_binned - 1

def TrackClassification(x,experiment):
	ev = np.ones(experiment.NumberOfEvents)
	ev[experiment.Sample==2] = x
	return experiment.Exp_wBinIt(ev) / experiment.weightOscBF_binned - 1

def IceAbsorption_ORCA(x, experiment):
	xx = x - 1
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return xx * d

def IceScattering_ORCA(x, experiment):
	xx = x - 1
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return xx * d

def OffSet_ORCA(x, experiment):
	xx = x
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return xx * d

def OptEffHeadon_ORCA(x, experiment):
	xx = x
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return xx * d

def OptEffLateral_ORCA(x, experiment):
	xx = x - 25
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return xx * d

def OptEffOverall_ORCA(x, experiment):
	xx = x - 1
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return xx * d

def CoinFraction_ORCA(x, experiment):
	xx = x
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return xx * d