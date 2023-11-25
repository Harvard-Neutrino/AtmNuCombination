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
