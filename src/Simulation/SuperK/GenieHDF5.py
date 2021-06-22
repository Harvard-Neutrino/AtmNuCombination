import h5py
import numpy as np
import ROOT as root
from applications import Azimuth
from math import log10


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

			# self.Flux()

	def GetLeptonPDG(self, ipnu):
		lep_pdg = ipnu
		lep_pdg[self.CC] = (abs(lep_pdg[self.CC]) - 1) * np.sign(lep_pdg[self.CC])
		self.PDGlep = lep_pdg

	# def RestHadrons(self, index):
	# 	return self.Phad[index]>0.001

	def Flux(self, tree): # To be implemented with IC's NuFlux
		CZ_points = 10
		E_points = 20
		Azi_step = 30*3.14159265/180.

		flux_nue = np.array([])
		flux_nueb = np.array([])
		flux_numu = np.array([])
		flux_numub = np.array([])

		E_per_Azi = 101
		Azi_per_CZ = 12

		for i,E in enumerate(self.Enu):
			E_index = int(E_points*(log10(E)+1))
			CZ_index = int(CZ_points*(1-self.Pznu[i]/E)) # decreasing order in cos zenith
			Azi_index = int(Azimuth(self.Pynu[i]/E, self.Pznu[i]/E) / Azi_step)

			tree_index0 = E_index + E_per_Azi*Azi_index + E_per_Azi*Azi_per_CZ*CZ_index
			tree_index1 = tree_index0 + 1

			tree.GetEntry(tree_index0)
			E0 = tree.logE
			flux0 = np.array([tree.flux[0],tree.flux[1],tree.flux[2],tree.flux[3]])
			tree.GetEntry(tree_index1)
			E1 = tree.logE
			flux1 = np.array([tree.flux[0],tree.flux[1],tree.flux[2],tree.flux[3]])
			
			flux_numu = np.append(flux_numu, np.interp(log10(E), [E0,E1], [flux0[0],flux1[0]]))
			flux_numub = np.append(flux_numub, np.interp(log10(E), [E0,E1], [flux0[1],flux1[1]]))
			flux_nue = np.append(flux_nue, np.interp(log10(E), [E0,E1], [flux0[2],flux1[2]]))
			flux_nueb = np.append(flux_nueb, np.interp(log10(E), [E0,E1], [flux0[3],flux1[3]]))


		self.Flux_numu = flux_numu
		self.Flux_numub = flux_numub
		self.Flux_nue = flux_nue
		self.Flux_nueb = flux_nueb


