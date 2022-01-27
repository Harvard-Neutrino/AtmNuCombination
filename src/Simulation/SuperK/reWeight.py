import numpy as np
import argparse
import h5py
import sys

def simMatrix(itype, nu, mode, oscw):
	# Build rate matrix for SuperKSim output
	cc = abs(mode)<30
	ccnue   = (cc==1)*(nu==12)
	ccnueb  = (cc==1)*(nu==-12)
	ccnumu  = (cc==1)*(abs(nu)==14)
	ccnutau = (cc==1)*(abs(nu)==16)
	nc      = (cc==0)

	matrix = np.zeros((16, 5))
	for t in range(16):
		a = (itype==t)*(ccnue)
		matrix[t,0] = np.sum(oscw[a])
		a = (itype==t)*(ccnueb)
		matrix[t,1] = np.sum(oscw[a])
		a = (itype==t)*(ccnumu)
		matrix[t,2] = np.sum(oscw[a])
		a = (itype==t)*(ccnutau)
		matrix[t,3] = np.sum(oscw[a])
		a = (itype==t)*(nc)
		matrix[t,4] = np.sum(oscw[a])

	norm = np.sum(matrix)
	if abs(np.sum(matrix)-norm)>0.0001:
		print('WARNING: Check the reweighting, there is a potential error here.')
	return matrix / norm

def allSKMatrix():
	matrix = np.array([[0.184519891405709,0.063822779733077,0.000514699836557,0.0,0.008492547303192],
	[0.023220765775251,0.000548067763639,0.003115332551214,2.88456717704981E-05,0.001932660008623],
	[0.001376049093042,0.000473016875733,0.000215007670788,0,0.012269771079626],
	[0.002902774017588,0.000920391761674,0.053736719008513,7.07993662826259E-05,0.013168682128569],
	[0.000200761563838,0,0.195140240050537,0,0.005420562223626],
	[0,0,0.016860015642391,1.7221670727672E-05,0.000344433414553],
	[0.002902936959159,0.001083185432522,0.00047660159031,0,0.038864693318894],
	[0.010450252684776,0.001514529374605,0.001682810416228,0.000555327437355,0.002625184249316],
	[0.030025373762898,0.020456848058238,0.000494923743344,0.000549915270383,0.003464466203411],
	[0.000193559546371,6.45198487902215E-05,0.0640036899999,0.00012903969758,0.00012903969758],
	[0.012782979374104,0.002340868754324,0.002685114159372,0.000917987746794,0.004222743635252],
	[0.010300108293475,0.005237343200072,0.000795300263715,0.000426746482969,0.002638069167444],
	[0.001552657702373,0.000230023363315,0.052502832676554,0.000287529204143,0.002932797882261],
	[0.012219740496746,0.001288627179657,0.015463526155883,0.00217733557942,0.013286190576462],
	[0.001200252684776,0.00045723911801,0.011845350900942,0.000142887224378,0.000642992509701],
	[0.000434662936558,0.000217331468279,0.070850058658966,0.000507106759318,0.000434662936558]])
	return matrix


def modeIndex(nu, mode):
	if abs(mode)<30:
		if nu==12:
			j=0
		elif nu==-12:
			j=1
		elif abs(nu)==14:
			j=2
		elif abs(nu)==16:
			j=3
		else:
			print('WARNING: out of range')
	else:
		j=4

	return j

