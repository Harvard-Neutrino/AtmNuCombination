import numpy as np
import h5py
from math import pi
import argparse


def SampleBins(itype):
	if itype==0 or itype==1:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
	elif itype==2:
		zentihBins = np.array([-1.0, 1.0])
		energyBins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
	elif itype>2 and itype<6:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
	elif itype==6:
		zentihBins = np.array([-1.0, 1.0])
		energyBins = np.array([0.1, 0.15, 0.25, 0.4, 0.63, 1.33])
	elif itype>6 and itype<9:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([1.3, 2.5, 5., 10., 100.])
	elif itype==9:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([1.3, 3.0, 100.])
	elif itype>9 and itype<12:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([1.3, 2.5, 5.0, 100.])
	elif itype==12:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([0.6, 1.3, 2.5, 5., 100.])
	elif itype==13:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([1.3, 2.5, 5.0, 10., 100.])
	elif itype==14:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([1.0, 20., 1.0e5])
	elif itype==15:
		zentihBins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
		energyBins = np.array([1.0, 10.0, 50., 1.0e3, 1.0e5])
	elif itype==16:
		zentihBins = np.array([-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
		energyBins = np.array([1.6, 10.0, 50., 1.0e5])
	elif itype>16:
		zentihBins = np.array([-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
		energyBins = np.array([1.6, 1.0e5])

	return zentihBins, energyBins

def SK_sample():
	s = {0:'SubGeV e-like 0dcy', 1:'SubGeV e-like 1dcy', 2:'Pi0-like 1ring', 3:'SubGeV mu-like 0dcy', 4:'SubGeV mu-like 1dcy', 5:'SubGeV mu-like 2dcy', 6:'Pi0-like 2ring', 7:'MultiGeV e-like nue', 8:'MultiGeV e-like nuebar', 9:'MultiGeV mu-like', 10:'MultiRing e-like nue', 11:'MultiRing e-like nuebar', 12:'MultiRing mu-like', 13:'MultiRing Other', 14:'PC Stoppping', 15:'PC Through-going', 16:'UpMu Stopping', 17:'UpMu non-Showering', 18:'UpMu Showering'}
	return s

def IC_sample():
	s = {0:'Tracks', 1:'No Tracks'}
	return s


parser = argparse.ArgumentParser()
parser.add_argument("unosc_hdf5filename", type=str, nargs='?', default='../Simulation/SuperK/data/testfcmc.hdf5')

parser.add_argument('theta23_min', type=float, nargs='?', default=0.3)
parser.add_argument('theta23_max', type=float, nargs='?', default=0.7)
parser.add_argument('theta23_points', type=int, nargs='?', default=10)

parser.add_argument('deltam32_min', type=float, nargs='?', default=0.0015)
parser.add_argument('deltam32_max', type=float, nargs='?', default=0.004)
parser.add_argument('deltam32_points', type=int, nargs='?', default=10)

parser.add_argument('theta13_min', type=float, nargs='?', default=0.018)
parser.add_argument('theta13_max', type=float, nargs='?', default=0.018)
parser.add_argument('theta13_points', type=int, nargs='?', default=1)

parser.add_argument('deltacp_min', type=float, nargs='?', default=0.0)
parser.add_argument('deltacp_max', type=float, nargs='?', default=2*pi)
parser.add_argument('deltacp_points', type=int, nargs='?', default=10)

args = parser.parse_args()
filename = args.unosc_hdf5filename

# Usually fixed parameters
th12 = 0.563942
dm21 = 7.65e-05

# Global fit free parameters
s2vTheta23  = np.linspace(args.theta23_min, args.theta23_max, num=args.theta23_points)
vTheta23  = np.arcsin(np.sqrt(s2vTheta23))
vDeltaM32 = np.geomspace(args.deltam32_min, args.deltam32_max, num=args.deltam32_points)
s2vTheta13  = np.linspace(args.theta13_min, args.theta13_max, num=args.theta13_points)
vTheta13  = np.arcsin(np.sqrt(s2vTheta13))
vDeltaCP  = np.linspace(args.deltacp_min, args.deltacp_max, num=args.deltacp_points)
vNMO      = np.array(['NO','IO'])
# print(s2vTheta23)
# print(s2vTheta23)
# print(vDeltaCP)
# print(vDeltaM32)

# Open reco file(s)
with h5py.File(filename, 'r') as hf:
	evis = np.array(hf['evis'][()])
	recocz = np.array(hf['recodirZ'][()])
	truecz = np.array(hf['dirnuZ'][()])
	mode = np.array(hf['mode'][()])
	ipnu = np.array(hf['ipnu'][()])
	pnu = np.array(hf['pnu'][()])
	itype = np.array(hf['itype'][()])
	dummy_osc = np.array(hf['oscw'][()])

# Experiment event samples
samples = SK_sample()

# H = 

with h5py.File('data/OscillatedBinnedData.hdf5','w') as out_hf:
	out_hf.create_dataset('Theta23', data=vTheta23, compression='gzip')
	out_hf.create_dataset('sin^2(Theta23)', data=s2vTheta23, compression='gzip')
	out_hf.create_dataset('Theta13', data=vTheta13, compression='gzip')
	out_hf.create_dataset('sin^2(Theta13)', data=s2vTheta13, compression='gzip')
	out_hf.create_dataset('DeltaM32(eV^2)', data=vDeltaM32, compression='gzip')
	out_hf.create_dataset('DeltaCP', data=vDeltaCP, compression='gzip')
	for th23 in s2vTheta23:
		for dm32 in vDeltaM32:
			for th13 in s2vTheta13:
				for dcp in vDeltaCP:
					for mo in vNMO:
						group = out_hf.create_group('/sin^2(Theta23)='+str(round(th23,2))+'/DeltaM32(eV^2)='+str(round(dm32,4))+'/sin^2(Theta13)='+str(round(th13,3))+'/DeltaCP='+str(round(dcp,2))+'/'+mo)
						
						for i in range(19):
							zbin, ebin = SampleBins(i)
							H, zbin, ebin = np.histogram2d(recocz[itype==i], evis[itype==i], bins=(zbin,ebin), weights=dummy_osc[itype==i])
							# print(H)
							group.create_dataset(samples[i]+' Events', data=H, compression='gzip')
							group.create_dataset(samples[i]+' Zenith bins', data=zbin, compression='gzip')
							group.create_dataset(samples[i]+' Energy bins', data=ebin, compression='gzip')
