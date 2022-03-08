import numpy as np
from Systematics.SKDetector import SKEnergyScale

####################
# Super-Kamiokande #
####################

def Diff_SKEnergyScale(x, experiment):
	h0 = x*(1+1e-6)
	h1 = x*(1-1e-6)
	w0 = SKEnergyScale(h0,experiment)
	w1 = SKEnergyScale(h1,experiment)
	dw = (w0 - w1) / (h0 - h1)
	return dw

def Diff_FCPCSeparation(x,experiment):
	fcpc = np.zeros(experiment.NumberOfEvents)
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		fcpc[experiment.Sample<16] = 1
		wFC = np.sum(experiment.Weight[experiment.Sample<16])
		wPC = np.sum(experiment.Weight[np.logical_or(experiment.Sample==16, experiment.Sample==17)])
		y = (-wFC) / wPC
		fcpc[np.logical_or(experiment.Sample==16,experiment.Sample==17)] = y
	else:
		fcpc[experiment.Sample<14] = 1
		wFC = np.sum(experiment.Weight[experiment.Sample<14])
		wPC = np.sum(experiment.Weight[np.logical_or(experiment.Sample==14,experiment.Sample==15)])
		y = (-wFC) / wPC
		fcpc[np.logical_or(experiment.Sample==14,experiment.Sample==15)] = y
	return fcpc

def Diff_FCReduction(x,experiment):
	fc = np.zeros(experiment.NumberOfEvents)
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		fc[experiment.Sample<16] = 1
	else:
		fc[experiment.Sample<14] = 1
	return fc

def Diff_FiducialVolume(x,experiment):
	return np.ones(experiment.NumberOfEvents)

def Diff_PCReduction(x,experiment):
	pc = np.zeros(experiment.NumberOfEvents)
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		pc[np.logical_or(experiment.Sample==16 , experiment.Sample==17)] = 1
	else:
		pc[np.logical_or(experiment.Sample==14, experiment.Sample==15)] = 1
	return pc

def Diff_SubGeV2ringPi0(x,experiment):
	pi02r = np.zeros(experiment.NumberOfEvents)
	pi02r[experiment.Sample==6] = 1
	return pi02r

def Diff_SubGeV1ringPi0(x,experiment):
	pi01r = np.zeros(experiment.NumberOfEvents)
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		pi01r[experiment.Sample==3] = 1
	else:
		pi01r[experiment.Sample==2] = 1
	return pi01r

def Diff_MultiRing_NuNuBarSeparation(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		nu = 12
		nub = 13
	else:
		nu = 10
		nub = 11
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = np.sum(experiment.Weight[experiment.Sample==nu])
	n1 = np.sum(experiment.Weight[experiment.Sample==nub])
	r = n0 / n1
	mr[experiment.Sample==nu] = 1
	mr[experiment.Sample==nub] = -r
	return mr

def Diff_MultiRing_EMuSeparation(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		e0 = 12
		e1 = 13
		e2 = 15
		mu = 14
	else:
		e0 = 10
		e1 = 11
		e2 = 13
		mu = 12
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = np.sum(experiment.Weight[experiment.Sample==e0]) + np.sum(experiment.Weight[experiment.Sample==e1]) + np.sum(experiment.Weight[experiment.Sample==e2])
	n1 = np.sum(experiment.Weight[experiment.Sample==mu])
	r = n0 / n1
	mr[experiment.Sample==e0] = 1
	mr[experiment.Sample==e1] = 1
	mr[experiment.Sample==e2] = 1
	mr[experiment.Sample==mu] = -r
	return mr


def Diff_MultiRing_EOtherSeparation(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		e0 = 12
		e1 = 13
		o0 = 15
	else:
		e0 = 10
		e1 = 11
		o0 = 13
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = np.sum(experiment.Weight[experiment.Sample==e0]) + np.sum(experiment.Weight[experiment.Sample==e1])
	n1 = np.sum(experiment.Weight[experiment.Sample==o0])
	r = n0 / n1
	mr[experiment.Sample==e0] = 1
	mr[experiment.Sample==e1] = 1
	mr[experiment.Sample==o0] = -r
	return mr

def Diff_PC_StopThruSeparation(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		pcs = 16
		pct = 17
	else:
		pcs = 14
		pct = 15
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = np.sum(experiment.Weight[experiment.Sample==pcs])
	n1 = np.sum(experiment.Weight[experiment.Sample==pct])
	r = n0 / n1
	mr[experiment.Sample==pcs] = 1
	mr[experiment.Sample==pct] = -r
	return mr

def Diff_Pi0_RingSeparation(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		r1  = 3
		r2  = 6
	else:
		r1  = 2
		r2  = 6
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = np.sum(experiment.Weight[experiment.Sample==r1 ])
	n1 = np.sum(experiment.Weight[experiment.Sample==r2 ])
	r = n0 / n1
	mr[experiment.Sample==r1 ] = 1
	mr[experiment.Sample==r2 ] = -r
	return mr

def Diff_E_RingSeparation(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		r1  = [0,1,2,7,8,9]
		r2  = [12,13,14]
	else:
		r1  = [0,1,2,7,8,9]
		r2  = [12,13,14]
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = 0
	n1 = 0
	for sample in r1 :
		n0 += np.sum(experiment.Weight[experiment.Sample==sample])
	for sample in r2 :
		n1 += np.sum(experiment.Weight[experiment.Sample==sample])
	r = n0 / n1
	for sample in r1 :
		mr[experiment.Sample==sample] = 1
	for sample in r2 :
		mr[experiment.Sample==sample] = -r
	return mr

def Diff_Mu_RingSeparation(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		r1  = [4,5,10,11]
		r2  = [14]
	else:
		r1  = [3,4,5,9]
		r2  = [12]
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = 0
	n1 = 0
	for sample in r1 :
		n0 += np.sum(experiment.Weight[experiment.Sample==sample])
	for sample in r2 :
		n1 += np.sum(experiment.Weight[experiment.Sample==sample])
	r = n0 / n1
	for sample in r1 :
		mr[experiment.Sample==sample] = 1
	for sample in r2 :
		mr[experiment.Sample==sample] = -r
	return mr

def Diff_SingleRing_PID(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		e = [0,1,2,3,7,8,9]
		mu = [4,5,10,11]
	else:
		e = [0,1,2,7,8]
		mu = [3,4,5,9]
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = 0
	n1 = 0
	for sample in e:
		n0 += np.sum(experiment.Weight[experiment.Sample==sample])
	for sample in mu:
		n1 += np.sum(experiment.Weight[experiment.Sample==sample])
	r = n0 / n1
	for sample in e:
		mr[experiment.Sample==sample] = 1
	for sample in mu:
		mr[experiment.Sample==sample] = -r
	return mr

def Diff_MultiRing_PID(x,experiment):
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		e = [12,13,15]
		mu = [14]
	else:
		e = [10,11,13]
		mu = [12]
	mr = np.zeros(experiment.NumberOfEvents)
	n0 = 0
	n1 = 0
	for sample in e:
		n0 += np.sum(experiment.Weight[experiment.Sample==sample])
	for sample in mu:
		n1 += np.sum(experiment.Weight[experiment.Sample==sample])
	r = n0 / n1
	for sample in e:
		mr[experiment.Sample==sample] = 1
	for sample in mu:
		mr[experiment.Sample==sample] = -r
	return mr

def Diff_NeutronTagging(x,experiment):
	nn = np.zeros(experiment.NumberOfEvents)
	if experiment.Experiment == 'SuperK-Gd' or experiment.Experiment == 'SKIV' or experiment.Experiment == 'SuperK_Htag' or experiment.Experiment == 'SuperK_Gdtag':
		n0 = np.sum(experiment.Neutron==0)
		n1 = np.sum(experiment.Neutron>0)
		r = n0 / n1
		nn[experiment.Neutron==0] = 1
		nn[experiment.Neutron>0] = -r
		return nn
	else:
		return 0

def Diff_DecayETagging(x,experiment):
	mue = np.zeros(experiment.NumberOfEvents)
	n0 = np.sum(experiment.DecayE==0)
	n1 = np.sum(experiment.DecayE==1)
	n2 = np.sum(experiment.DecayE>1)
	N = n0 + n1 + n2
	r0 = n0 / N
	r1 = n1 / N
	r2 = n2 / N
	rx1 = r1 - 2*r2
	rx2 = 2*x*r2 - 2*r2
	rx0 = - rx1 - rx2
	mue[experiment.DecayE==0] = rx0 / r0
	mue[experiment.DecayE==1] = rx1 / r1
	mue[experiment.DecayE>1] = rx2 / r2
	return mue
