import numpy as np
import argparse
import h5py
import sys
from reWeight import *

parser = argparse.ArgumentParser()
parser.add_argument('fname', type=str, nargs='?', default='NULL')
args = parser.parse_args()
infile = args.fname

if infile == 'NULL':
    sys.exit('Please introduce your HDF5 file')

with h5py.File(infile, 'r') as hf:
	mode = np.array(hf['mode'][()])
	ipnu = np.array(hf['ipnu'][()])
	weightOsc_SKpaper = np.array(hf['weightOsc_SKpaper'][()])
	itype = np.array(hf['itype'][()])

W = simMatrix(itype, ipnu, mode, weightOsc_SKpaper) # Rate matrix from this simulation
W0= allSKMatrix() # Rate matrix from SK's paper
weightReco = np.zeros(np.size(weightOsc_SKpaper))
for i,t in enumerate(itype):
	ti = int(t)
	j = modeIndex(ipnu[i], mode[i])
	if ti>-1 and ti<16:
		if W[ti][j]>0:
			weightReco[i] = W0[ti][j] / W[ti][j]
		else:
			print('Simulation weight is zero!')
			weightReco[i] = 0
	else:
		weightReco[i] = 0

with h5py.File(infile, 'r+') as hf:
	del hf['weightReco']
	hf.create_dataset('weightReco', data=weightReco, compression='gzip')
