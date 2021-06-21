import h5py
import numpy as np

class GenieSimulation:
	def __init__(self, filename):
		# allGenieKeys = np.array(['A', 'Ef', 'Ei', 'El', 'En', 'Ev', 'EvRF', 'Q2', 'Q2s', 'W', 'Ws', 'Z', 'calresp0', 
		# 	'cc', 'charm', 'coh', 'cthf', 'cthl', 'dfr', 'dis', 'em', 'fspl', 'hitnuc', 'hitqrk', 'iev', 'imd', 'imdanh', 
		# 	'mec', 'nc', 'neu', 'neut_code', 'nf', 'nfem', 'nfk0', 'nfkm', 'nfkp', 'nfn', 'nfother', 'nfp', 'nfpi0', 
		# 	'nfpim', 'nfpip', 'ni', 'niem', 'nik0', 'nikm', 'nikp', 'nin', 'niother', 'nip', 'nipi0', 'nipim', 'nipip', 
		# 	'nuance_code', 'nuel', 'pdgf', 'pdgi', 'pf', 'pl', 'pxf', 'pxi', 'pxl', 'pxn', 'pxv', 'pyf', 'pyi', 'pyl', 
		# 	'pyn', 'pyv', 'pzf', 'pzi', 'pzl', 'pzn', 'pzv', 'qel', 'res', 'resc', 'resid', 'sea', 'singlek', 'sumKEf', 
		# 	't', 'tgt', 'ts', 'vtxt', 'vtxx', 'vtxy', 'vtxz', 'wght', 'x', 'xs', 'y', 'ys'])
		usedGenieKeys = np.array(['Ef', 'El', 'Ev', 'cc', 'nc', 'neu', 'neut_code', 'nf', 'pdgf', 
			'pf', 'pl', 'pxf', 'pxl', 'pxv', 'pyf', 'pyl', 'pyv', 'pzf', 'pzl', 'pzv'])

		with h5py.File(filename, 'r') as hf:
			self.Enu = np.array(hf['Ev'][()])
			self.Ipnu = np.array(hf['neu'][()])
			self.CC = np.array(hf['cc'][()])
			self.NC = np.array(hf['nc'][()])
			self.Mode = np.array(hf['neut_code'][()])
			self.Pxnu = np.array(np.array(hf['pxv'][()]))
			self.Pynu = np.array(np.array(hf['pyv'][()]))
			self.Pznu = np.array(np.array(hf['pzv'][()]))
			
			self.Elep = np.array(hf['El'][()])
			self.GetLeptonPDG(np.array(hf['neu'][()]))
			self.Pxlep = np.array(np.array(hf['pxl'][()]))
			self.Pylep = np.array(np.array(hf['pyl'][()]))
			self.Pzlep = np.array(np.array(hf['pzl'][()]))
			self.Plep = np.array(np.array(hf['pl'][()]))

			self.Nhad = np.array(hf['nf'][()])
			self.Ehad = np.array(hf['Ef'][()])
			self.Phad = np.array(hf['pf'][()])
			self.PDGhad = np.array(hf['pdgf'][()])
			self.Pxhad = np.array(np.array(hf['pxf'][()]))
			self.Pyhad = np.array(np.array(hf['pyf'][()]))
			self.Pzhad = np.array(np.array(hf['pzf'][()]))

			self.Flux()

	def GetLeptonPDG(self, ipnu):
		lep_pdg = ipnu
		lep_pdg[self.CC] = (abs(lep_pdg[self.CC]) - 1) * np.sign(lep_pdg[self.CC])
		self.PDGlep = lep_pdg

	# def RestHadrons(self, index):
	# 	return self.Phad[index]>0.001

	def Flux(self):
		flux = np.ones(4)
		self.Flux = flux
