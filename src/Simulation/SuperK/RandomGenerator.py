import h5py
import numpy as np
from numpy.random import choice

class RecoDists:
	def __init__(self, filename='defaultRecoDistributions.hdf5'):
		self.f = h5py.File(filename, 'r')

	def Random(self, label):
		w = np.array(self.f[label]['weights'])
		v = np.array(self.f[label]['values'])
		value = choice(v, p=w)
		return value

	def Random_and_Weight(self, label):
		w = np.array(self.f[label]['weights'])
		v = np.array(self.f[label]['values'])
		value = choice(v, p=w)
		return value, w[np.where(v==value)]



