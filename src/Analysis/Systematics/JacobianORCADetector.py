import numpy as np

####################
####### ORCA #######
####################

def Diff_ORCAEnergyScale(x, experiment):
	if np.abs(x-1)<1e-3: 
		return 0
	h0 = x+1e-3
	h1 = x-1e-3
	w0 = experiment.Exp_wBinIt(1,shift_E=h0)
	w1 = experiment.Exp_wBinIt(1,shift_E=h1)
	dw = (w0 - w1) / (h0 - h1) / experiment.weightOscBF_binned
	return dw
	
def Diff_ShowerClassification(x,experiment):
	ev = np.zeros(experiment.NumberOfEvents)
	ev[experiment.Sample==0] = 1
	return experiment.Exp_wBinIt(ev) / experiment.weightOscBF_binned

def Diff_IntermediateClassification(x,experiment):
	ev = np.zeros(experiment.NumberOfEvents)
	ev[experiment.Sample==1] = 1
	return experiment.Exp_wBinIt(ev) / experiment.weightOscBF_binned

def Diff_TrackClassification_ORCA(x, experiment):
	ev = np.zeros(experiment.NumberOfEvents)
	ev[experiment.Sample==2] = 1
	return experiment.Exp_wBinIt(ev) / experiment.weightOscBF_binned

def Diff_IceAbsorption_ORCA(x, experiment):
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return d

def Diff_IceScattering_ORCA(x, experiment):
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return d

def Diff_OffSet_ORCA(x, experiment):
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return d

def Diff_OptEffHeadon_ORCA(x, experiment):
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return d

def Diff_OptEffLateral_ORCA(x, experiment):
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return d

def Diff_OptEffOverall_ORCA(x, experiment):
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return d

def Diff_CoinFraction_ORCA(x, experiment):
	d = experiment.ExpFracNuECC * experiment.ice_absorption['nueCC'] + \
	experiment.ExpFracNuMuCC * experiment.ice_absorption['numuCC'] + \
	experiment.ExpFracNuTauCC * experiment.ice_absorption['nutauCC'] + \
	experiment.ExpFracNC * experiment.ice_absorption['NC']
	return d
