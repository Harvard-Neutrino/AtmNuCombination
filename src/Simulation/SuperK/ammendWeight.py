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
	itype = np.array(hf['otype'][()])

# First tune to match table
W = simMatrix(itype, ipnu, mode, weightOsc_SKpaper) # Rate matrix from this simulation
W0= allSKMatrix() # Rate matrix from SK's paper

weightReco = np.ones(np.size(weightOsc_SKpaper))


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

# Transfer events between samples if weightReco found too large
SampleGroups = {0:[0,1], 1:[2,6], 2:[3,4,5], 3:[7,8], 4:[9,12], 5:[10,11,13], 6:[14,15]}
NuMode = {0:(ipnu==12) * (np.abs(mode)<30), 1:(ipnu==-12) * (np.abs(mode)<30), 2:(np.abs(ipnu)==14) * (np.abs(mode)<30), 3:(np.abs(ipnu)==16) * (np.abs(mode)<30), 4:(np.abs(mode)>30)}
t_itype = itype
for gr in SampleGroups:
	for m in NuMode:
		p = []
		cond = 0
		for sample in SampleGroups[gr]:
			cond = NuMode[m] * (itype==sample)
			index = np.where(cond)
			index = np.ravel(index)
			for a in index:
				s=a
				break
			wOr= weightReco[s]
			if wOr>0:
				prob = 1-wOr
				# prob = 1-1/wOr
				if prob>1:
					prob = 1
				elif prob<-1:
					prob = -1
			if wOr==0:
				prob = 1

			print(f'Sample {sample}, mode {m}, weight {wOr}, prob {prob}')
			p.append(prob) # Transfer probability

		sample = SampleGroups[gr]
		if len(p)==2:
			if p[0]>0 and p[1]<0:
				print('check above')
				cond = (NuMode[m]) * (itype==sample[0])
				index = np.where(cond)
				index = index[0]
				index2transfer = np.random.choice(index, int(p[0]*index.size))
				t_itype[index2transfer] = sample[1]
			elif p[1]>0 and p[0]<0:
				print('check above')
				cond = (NuMode[m]) * (itype==sample[1])
				index = np.where(cond)
				index = index[0]
				index2transfer = np.random.choice(index, int(p[1]*index.size))
				t_itype[index2transfer] = sample[0]

		elif len(p)==3:
			if p[0]>0:
				if p[1]<0 and p[2]<0: # 0 --> 1,2
					print('check above')
					cond = (NuMode[m]) * (itype==sample[0])
					index = np.where(cond)
					index = index[0]
					frac = p[1] / (p[1]+p[2])
					index2transfer1 = np.random.choice(index, int(frac*p[0]*index.size))
					index2transfer2 = np.random.choice(index, int((1-frac)*p[0]*index.size))
					t_itype[index2transfer1] = sample[1]
					t_itype[index2transfer2] = sample[2]
				elif p[1]>0 and p[2]<0: # 0,1 --> 2
					print('check above')
					cond = (NuMode[m]) * (itype==sample[0])
					index = np.where(cond)
					index = index[0]
					index2transfer1 = np.random.choice(index, int(p[0]*index.size))
					cond = (NuMode[m]) * (itype==sample[1])
					index = np.where(cond)
					index = index[0]
					index2transfer2 = np.random.choice(index, int(p[1]*index.size))
					t_itype[index2transfer1] = sample[2]
					t_itype[index2transfer2] = sample[2]
				elif p[1]<0 and p[2]>0: # 0,2 --> 1
					print('check above')
					cond = (NuMode[m]) * (itype==sample[0])
					index = np.where(cond)
					index = index[0]
					index2transfer1 = np.random.choice(index, int(p[0]*index.size))
					cond = (NuMode[m]) * (itype==sample[2])
					index = np.where(cond)
					index = index[0]
					index2transfer2 = np.random.choice(index, int(p[2]*index.size))
					t_itype[index2transfer1] = sample[1]
					t_itype[index2transfer2] = sample[1]
			elif p[1]>0:
				if p[0]<0 and p[2]<0: # 1 --> 0,2
					print('check above')
					cond = (NuMode[m]) * (itype==sample[1])
					index = np.where(cond)
					index = index[0]
					frac = p[0] / (p[0]+p[2])
					index2transfer1 = np.random.choice(index, int(frac*p[1]*index.size))
					index2transfer2 = np.random.choice(index, int((1-frac)*p[1]*index.size))
					t_itype[index2transfer1] = sample[0]
					t_itype[index2transfer2] = sample[2]
				elif p[0]<0 and p[2]>0: # 1,2 --> 0
					print('check above')
					cond = (NuMode[m]) * (itype==sample[1])
					index = np.where(cond)
					index = index[0]
					index2transfer1 = np.random.choice(index, int(p[1]*index.size))
					cond = (NuMode[m]) * (itype==sample[2])
					index = np.where(cond)
					index = index[0]
					index2transfer2 = np.random.choice(index, int(p[2]*index.size))
					t_itype[index2transfer1] = sample[0]
					t_itype[index2transfer2] = sample[0]
			elif p[2]>0:
				if p[0]<0 and p[2]<0: # 2 --> 0,1
					print('check above')
					cond = (NuMode[m]) * (itype==sample[2])
					index = np.where(cond)
					index = index[0]
					frac = p[0] / (p[0]+p[1])
					index2transfer1 = np.random.choice(index, int(frac*p[2]*index.size))
					index2transfer2 = np.random.choice(index, int((1-frac)*p[2]*index.size))
					t_itype[index2transfer1] = sample[0]
					t_itype[index2transfer2] = sample[1]



# Second tune to match table

W = simMatrix(t_itype, ipnu, mode, weightOsc_SKpaper) # Rate matrix from this simulation
W0= allSKMatrix() # Rate matrix from SK's paper
dweightReco = np.zeros(np.size(weightOsc_SKpaper))


for i,t in enumerate(t_itype):
	ti = int(t)
	j = modeIndex(ipnu[i], mode[i])
	if ti>-1 and ti<16:
		if W[ti][j]>0:
			dweightReco[i] = W0[ti][j] / W[ti][j]
		else:
			print('Simulation weight is zero!')
			dweightReco[i] = 0
	else:
		dweightReco[i] = 0
'''
# Check effects
print('=========================================================')
for gr in SampleGroups:
	for m in NuMode:
		p = []
		cond = 0
		for sample in SampleGroups[gr]:
			cond = NuMode[m] * (t_itype==sample)
			index = np.where(cond)
			index = np.ravel(index)
			for a in index:
				s=a
				break
			wOr= dweightReco[s]
			if wOr>0:
				prob = 1-1/wOr
				if prob>1:
					prob = 1
				elif prob<-1:
					prob = -1
			if wOr==0:
				prob = 1

			print(f'Sample {sample}, mode {m}, weight {wOr}, prob {prob}')
			p.append(prob) # Transfer probability
'''


with h5py.File(infile, 'r+') as hf:
	del hf['weightReco']
	del hf['otype']
	del hf['itype']
	hf.create_dataset('weightReco', data=dweightReco, compression='gzip')
	hf.create_dataset('otype', data=t_itype, compression='gzip')
	hf.create_dataset('itype', data=t_itype, compression='gzip')
