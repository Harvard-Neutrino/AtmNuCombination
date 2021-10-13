import numpy as np
import h5py
import argparse

# The output of this code is analogous to the input hdf5 file, but with X2 values instead of number of entries

def poissonChi2(truth, test): # input: 2 2D histograms. output: 1 2d histogram
	# Compute statistical chi2 assuming Poisson-distributed data
	A = np.divide(test,truth)
	logA = np.log(A)
	OlogA = np.multiply(test,logA)
	chi2 = 2 * (truth - test + OlogA)
	return chi2

def SK_sample():
	s = {0:'SubGeV e-like 0dcy', 1:'SubGeV e-like 1dcy', 2:'Pi0-like 1ring', 3:'SubGeV mu-like 0dcy', 4:'SubGeV mu-like 1dcy', 5:'SubGeV mu-like 2dcy', 6:'Pi0-like 2ring', 7:'MultiGeV e-like nue', 8:'MultiGeV e-like nuebar', 9:'MultiGeV mu-like', 10:'MultiRing e-like nue', 11:'MultiRing e-like nuebar', 12:'MultiRing mu-like', 13:'MultiRing Other', 14:'PC Stoppping', 15:'PC Through-going', 16:'UpMu Stopping', 17:'UpMu non-Showering', 18:'UpMu Showering'}
	return s

# def IC_sample():
# 	s = {0:'Tracks', 1:'No Tracks'}
# 	return s

parser = argparse.ArgumentParser()
parser.add_argument('hdf5_oschisto_file', type=str, nargs='?', default='data/OscillatedBinnedData.hdf5')
parser.add_argument('hdf5_chi2histo_file', type=str, nargs='?', default='data/OscillatedBinnedChi2.hdf5')
parser.add_argument('true_s2th23', type=float, nargs='?', default='0.5')
parser.add_argument('true_dm32', type=float, nargs='?', default='0.0025')
parser.add_argument('true_dcp', type=float, nargs='?', default='4.1')
parser.add_argument('true_s2th13', type=float, nargs='?', default='0.018')
parser.add_argument('true_ordering', type=str, nargs='?', default='NO')

args = parser.parse_args()
infile = args.hdf5_oschisto_file
outfile = args.hdf5_chi2histo_file

# Will calculate the chi2 w.r.t this point
true_th23 = args.true_s2th23
true_th13 = args.true_s2th13
true_dm32 = args.true_dm32
true_dcp  = args.true_dcp
true_nmo  = args.true_ordering

samples = SK_sample()

with h5py.File(outfile, 'w') as ohf:
	with h5py.File(infile, 'r') as hf:
		s2vTheta23 = np.array(hf['sin^2(Theta23)'][()])
		s2vTheta13 = np.array(hf['sin^2(Theta13)'][()])
		vDeltaCP = np.array(hf['DeltaCP'][()])
		vDeltaM32 = np.array(hf['DeltaM32(eV^2)'][()])
		vNMO = np.array(['NO','IO'])
		ohf.create_dataset('sin^2(Theta23)', data=s2vTheta23, compression='gzip')
		ohf.create_dataset('sin^2(Theta13)', data=s2vTheta13, compression='gzip')
		ohf.create_dataset('DeltaM32(eV^2)', data=vDeltaM32, compression='gzip')
		ohf.create_dataset('DeltaCP', data=vDeltaCP, compression='gzip')

		#Pick closest values to inputs
		th23 = s2vTheta23[np.abs(s2vTheta23-true_th23).argmin()]
		th13 = s2vTheta13[np.abs(s2vTheta13-true_th13).argmin()]
		dcp  = vDeltaCP[np.abs(vDeltaCP-true_dcp).argmin()]
		dm32 = vDeltaM32[np.abs(vDeltaM32-true_dm32).argmin()]
		print('--------------------------')
		print('Assuming true point is:')
		print('--------------------------')
		print('sin^2(theta_23) = '+str(round(th23,2)))
		print('sin^2(theta_13) = '+str(round(th13,3)))
		print('Delta m^2_32 (eV^2) = '+str(round(dm32,4)))
		print('Delta_CP = '+str(round(dcp,2)))
		print('Ordering = '+true_nmo)
		print('===========================')

		# Build path to truth histograms
		truth_path = '/sin^2(Theta23)='+str(round(th23,2))+'/DeltaM32(eV^2)='+str(round(dm32,4))+'/sin^2(Theta13)='+str(round(th13,3))+'/DeltaCP='+str(round(dcp,2))+'/'+true_nmo
		truth_location = hf[truth_path]
		# Go through the rest of oscillation points
		for test_th23 in s2vTheta23:
			for test_dm32 in vDeltaM32:
				for test_th13 in s2vTheta13:
					for test_dcp in vDeltaCP:
						for test_mo in vNMO:
							test_path = '/sin^2(Theta23)='+str(round(test_th23,2))+'/DeltaM32(eV^2)='+str(round(test_dm32,4))+'/sin^2(Theta13)='+str(round(test_th13,3))+'/DeltaCP='+str(round(test_dcp,2))+'/'+test_mo
							test_location = hf[test_path]
							group = ohf.create_group(test_path)
							for i in samples:
								truth_H = np.array(truth_location[samples[i]+' Events'][()])
								cz = np.array(truth_location[samples[i]+' Zenith bins'][()])
								ev = np.array(truth_location[samples[i]+' Energy bins'][()])
								test_H = np.array(test_location[samples[i]+' Events'][()])
								x2 = poissonChi2(truth_H, test_H)
								group.create_dataset(samples[i]+' chi2 (poisson stats)', data=x2, compression='gzip')
								group.create_dataset(samples[i]+' Zenith bins', data=cz, compression='gzip')
								group.create_dataset(samples[i]+' Energy bins', data=ev, compression='gzip')



			
