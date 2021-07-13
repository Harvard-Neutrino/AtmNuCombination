import h5py
import numpy as np
import nuflux
import nuSQuIDS as nsq
import nuSQUIDSTools
import math

import sys

def FluxFactor(i, flavor, nue, nueb, numu, numub):
    j = int(abs(flavor) / 2) % 6
    factor = 1.0
    if i==j or i==1 and j==2:
        factor = 1.0
    elif i==0 and j>=1:
        if flavor>0:
            factor = nue / numu
        else:
            factor = nueb / numub
    elif i==1 and j==0:
        if flavor>0:
            factor = numu / nue
        else:
            factor = numub / nueb
    
    return factor


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
			self.Dirxnu = self.Pxnu / self.Enu
			self.Dirynu = self.Pynu / self.Enu
			self.Dirznu = self.Pznu / self.Enu
			self.Cz = self.Dirznu
			self.Azimuth()
			
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
		self.PointOsc()
		# self.MRPCUMflag()

	def MRPCUMflag(self):
		pass

	def Azimuth(self):
		self.Azi = np.arcsin(self.Dirynu/(np.sqrt(1-self.Cz**2))) + math.pi

	def GetLeptonPDG(self, ipnu):
		lep_pdg = ipnu
		lep_pdg[self.CC] = (abs(lep_pdg[self.CC]) - 1) * np.sign(lep_pdg[self.CC])
		self.PDGlep = lep_pdg

	def Flux(self): # To be implemented with IC's NuFlux
		# Try other fluxes, may need to tune nuflux to extend energy regions and include fluxes at SK location
		# flux = nuflux.makeFlux('honda2006')
		flux = nuflux.makeFlux('IPhonda2014_sk_solmin')
		numu = nuflux.NuMu
		numub = nuflux.NuMuBar
		nue = nuflux.NuE
		nueb = nuflux.NuEBar

		flux_nue = np.array([])
		flux_nueb = np.array([])
		flux_numu = np.array([])
		flux_numub = np.array([])

		for i,E in enumerate(self.Enu):
			# flux_numu  = np.append(flux_numu, flux.getFlux(numu, E, self.Azi[i], self.Cz[i]))
			# flux_numub = np.append(flux_numub, flux.getFlux(numub, E, self.Azi[i], self.Cz[i]))
			# flux_nue   = np.append(flux_nue, flux.getFlux(nue, E, self.Azi[i], self.Cz[i]))
			# flux_nueb  = np.append(flux_nueb, flux.getFlux(nueb, E, self.Azi[i], self.Cz[i]))
			flux_numu  = np.append(flux_numu, flux.getFlux(numu, E, self.Cz[i]))
			flux_numub = np.append(flux_numub, flux.getFlux(numub, E, self.Cz[i]))
			flux_nue   = np.append(flux_nue, flux.getFlux(nue, E, self.Cz[i]))
			flux_nueb  = np.append(flux_nueb, flux.getFlux(nueb, E, self.Cz[i]))

		self.Flux_numu = flux_numu
		self.Flux_numub = flux_numub
		self.Flux_nue = flux_nue
		self.Flux_nueb = flux_nueb

	def PointOsc(self):
		self.oscw = np.array([])
		units = nsq.Const()

		for k,(nu,E,cz,mod) in enumerate(zip(self.Ipnu, self.Enu, self.Cz, self.Mode)):
		# Get P_{x->ipnu} probabilities
			weight = 0.0
			if nu>0:
				nuSQ = nsq.nuSQUIDS(3,nsq.NeutrinoType.neutrino)
			elif nu<0:
				nuSQ = nsq.nuSQUIDS(3,nsq.NeutrinoType.antineutrino)
			else:
				print('What?! No identified neutrino flavour')
			nuSQ.Set_E(E*units.GeV)
			zenith = np.arccos(cz)
			nuSQ.Set_Track(nsq.EarthAtm().Track(zenith))
			nuSQ.Set_Body(nsq.EarthAtm())
			nuSQ.Set_rel_error(1.0e-4);
			nuSQ.Set_abs_error(1.0e-4);
			nuSQ.Set_MixingAngle(0,1,0.563942)
			nuSQ.Set_MixingAngle(0,2,math.asin(math.sqrt(0.018)))
			nuSQ.Set_MixingAngle(1,2,math.asin(math.sqrt(0.5)))
			nuSQ.Set_SquareMassDifference(1,7.65e-05)
			nuSQ.Set_SquareMassDifference(2,0.00247)
			nuSQ.Set_CPPhase(0,2,0)

			if abs(mod) < 30:
				for i in range(2):
					in_state = np.zeros(3)
					in_state[i] = 1
					Ffactor = FluxFactor(i, nu, self.Flux_nue[k], self.Flux_nueb[k], self.Flux_numu[k], self.Flux_numub[k])
					nuSQ.Set_initial_state(in_state,nsq.Basis.flavor)
					nuSQ.EvolveState()
					j = int(abs(nu) / 2) % 6
					prob = nuSQ.EvalFlavor(j)
					weight += prob*Ffactor
			else:
				weight = 1.0

			self.oscw = np.append(self.oscw,weight)
		print('Done with oscillations')


