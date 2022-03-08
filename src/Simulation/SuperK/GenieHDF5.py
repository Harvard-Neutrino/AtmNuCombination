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

		self.TopologySample()
		self.FluxWeight()
		self.Flux()
		self.PointOsc_SKTable()
		self.PointOsc_SKBest()

	def TopologySample(self):
		# SK topologies to choose from
		skTopology = np.array(['FC','PC-Stop','PC-Thru','UpMu-Stop','UpMu-Thru','UpMu-Show'])
		dummySample = np.array([])

		loge = np.zeros(60)
		line = np.zeros(60)
		fce  = np.zeros(60)
		fcm  = np.zeros(60)
		pcs  = np.zeros(60)
		pct  = np.zeros(60)
		ums  = np.zeros(60)
		umt  = np.zeros(60)
		umsh = np.zeros(60)
		# Acquiring digitized data
		with open('lib/SKTopologyFraction.dat') as f:
			lines = f.readlines()
			for i,l in enumerate(lines):
				loge[i], line[i], fce[i], fcm[i], pcs[i], pct[i], ums[i], umt[i], umsh[i] = l.split( )

		# CC electronic
		nue = fce + 0.116*pcs + 0.009*pct + 2*(0.011*ums + 0.003*umt + 0.001*umsh) # Factor x2 accounts for the upward-going cut of UpMus to be applied later
		fc_nue   = fce / nue
		pcs_nue  = 0.116*pcs / nue
		pct_nue  = 0.009*pct / nue
		ums_nue  = 2*0.011*ums / nue
		umt_nue  = 2*0.003*umt / nue
		umsh_nue = 2*0.001*umsh / nue
		# CC muonic
		numu = fcm + 0.829*pcs + 0.978*pct + 2*(0.986*ums + 0.996*umt + 0.998*umsh)
		numu[numu==0]=1.
		fc_numu   = fcm / numu
		pcs_numu  = 0.829*pcs / numu
		pct_numu  = 0.978*pct / numu
		ums_numu  = 2*0.986*ums / numu
		umt_numu  = 2*0.996*umt / numu
		umsh_numu = 2*0.998*umsh / numu
		# CC tauonic
		nut = 0.0057*(fce+fcm) + 0.01*pcs + 0.007*pct + 2*(0.0*ums + 0.0*umt + 0.0*umsh)
		nut[nut==0] = 1.
		fc_nut   = 0.0057*(fce+fcm) / nut
		pcs_nut  = 0.01*pcs / nut
		pct_nut  = 0.007*pct / nut
		ums_nut  = 2*0.0*ums / nut
		umt_nut  = 2*0.0*umt / nut
		umsh_nut = 2*0.0*umsh / nut
		# NC allic
		nc = 0.12*(fce+fcm) + 0.045*pcs + 0.006*pct + 2*(0.003*ums + 0.001*umt + 0.001*umsh)
		fc_nc   = 0.12*(fce+fcm) / nc
		pcs_nc  = 0.045*pcs / nc
		pct_nc  = 0.006*pct / nc
		ums_nc  = 2*0.003*ums / nc
		umt_nc  = 2*0.001*umt / nc
		umsh_nc = 2*0.001*umsh / nc

		# UpMu: Later cut on reconstructed direction
		for nu, E, cc in zip(self.Ipnu, self.Enu, self.CC):
			# Number of energy bin
			k = int((math.log10(E)+1)/0.1)
			if cc:
				if abs(nu)==12:
					probTopo = np.array([fc_nue[k],pcs_nue[k],pct_nue[k],ums_nue[k],umt_nue[k],umsh_nue[k]])
				elif abs(nu)==14:
					probTopo = np.array([fc_numu[k],pcs_numu[k],pct_numu[k],ums_numu[k],umt_numu[k],umsh_numu[k]])
				elif abs(nu)==16:
					probTopo = np.array([fc_nut[k],pcs_nut[k],pct_nut[k],ums_nut[k],umt_nut[k],umsh_nut[k]])
				else:
					print('WTF!?')
			else:
				probTopo = np.array([fc_nut[k],pcs_nut[k],pct_nut[k],ums_nut[k],umt_nut[k],umsh_nut[k]])

			if np.isnan(np.sum(probTopo)) or np.sum(probTopo)==0:
				sample = 'None'
			else:
				sample = np.random.choice(skTopology,1,p=probTopo)
			dummySample = np.append(dummySample, sample)
		self.TopologySample = dummySample


	def Azimuth(self):
		self.Azi = np.arcsin(self.Dirynu/(np.sqrt(1-self.Cz**2))) + 0.5*math.pi

	def GetLeptonPDG(self, ipnu):
		lep_pdg = ipnu
		lep_pdg[self.CC] = (abs(lep_pdg[self.CC]) - 1) * np.sign(lep_pdg[self.CC])
		self.PDGlep = lep_pdg

	def Flux(self):
		# flux = nuflux.makeFlux('IPhonda2014_sk_solmax')
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

	def FluxWeight(self): # Computes the inverse of the simulated flux used in the GENIE production for a given neutrino flavour
		flux = nuflux.makeFlux('IPhonda2014_sk_solmin')
		# flux = nuflux.makeFlux('IPhonda2014_sk_solmax')
		nus = {12:nuflux.NuE, -12:nuflux.NuEBar, 14:nuflux.NuMu, -14:nuflux.NuMuBar, 16:nuflux.NuMu, -16:nuflux.NuMuBar}
		flx_weight = np.array([])

		for i,E in enumerate(self.Enu):
			flx_weight  = np.append(flx_weight, 1./flux.getFlux(nus[self.Ipnu[i]], E, self.Cz[i]))

		self.FluxWeight = flx_weight

	def AtmInitialFlux(self, energies, zeniths, neutrino_flavors):
		AtmFlux = np.zeros((len(zeniths),len(energies),2,neutrino_flavors))
		flux = nuflux.makeFlux('IPhonda2014_sk_solmin')
		# flux = nuflux.makeFlux('IPhonda2014_sk_solmax')
		for ic,cz in enumerate(zeniths):
			for ie,E in enumerate(energies):
				AtmFlux[ic][ie][0][0] = flux.getFlux(nuflux.NuE,E,cz) # nue
				AtmFlux[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,E,cz) # nue bar
				AtmFlux[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,E,cz) # numu
				AtmFlux[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,E,cz) # numu bar
				AtmFlux[ic][ie][0][2] = 0.0 # nutau
				AtmFlux[ic][ie][1][2] = 0.0 # nutau bar
		return AtmFlux

		
	def PointOsc_SKBest(self):
		self.weightOsc_SKbest = np.array([])
		units = nsq.Const()

		for k,(nu,E,cz,mod) in enumerate(zip(self.Ipnu, self.Enu, self.Cz, self.Mode)):
		# Get P_{x->ipnu} probabilities
			weight = 0.0
			if nu>0:
				rho = 0
				# nuSQ = nsq.nuSQUIDS(3,nsq.NeutrinoType.neutrino)
			elif nu<0:
				rho = 1
				# nuSQ = nsq.nuSQUIDS(3,nsq.NeutrinoType.antineutrino)
			else:
				print('What?! No identified neutrino flavour')

			if cz>0:
				coszen = np.array([0.99999*cz,1.00001*cz])
			else:
				coszen = np.array([1.00001*cz,0.99999*cz])

			ener = np.array([0.99*E,1.01*E])
			AtmOsc = nsq.nuSQUIDSAtm(coszen,ener*units.GeV,3,nsq.NeutrinoType.both,False)
			AtmOsc.Set_rel_error(1.0e-5);
			AtmOsc.Set_abs_error(1.0e-5);
			AtmOsc.Set_MixingAngle(0,1,math.asin(math.sqrt(0.304)))
			AtmOsc.Set_MixingAngle(0,2,math.asin(math.sqrt(0.0219)))
			AtmOsc.Set_MixingAngle(1,2,math.asin(math.sqrt(0.588)))
			AtmOsc.Set_SquareMassDifference(1,7.53e-05)
			AtmOsc.Set_SquareMassDifference(2,0.0025)
			AtmOsc.Set_CPPhase(0,2,4.18)

			if abs(mod) < 30:
				for i in range(2):
					in_state = np.zeros((2,2,2,3))
					for ii in range(2):
						for ij in range(2):
							for ik in range(2):
								in_state[ii][ij][ik][i] = 1
					Ffactor = FluxFactor(i, nu, self.Flux_nue[k], self.Flux_nueb[k], self.Flux_numu[k], self.Flux_numub[k])
					AtmOsc.Set_initial_state(in_state,nsq.Basis.flavor)
					AtmOsc.EvolveState()
					j = int(abs(nu) / 2) % 6
					prob = 0
					for _ in range(20):
						prob += AtmOsc.EvalFlavor(j, cz, E*units.GeV, rho, True)
					prob /= 20
					weight += prob*Ffactor
			else:
				weight = 1.0

			self.weightOsc_SKbest = np.append(self.weightOsc_SKbest,weight)
		print('Done with oscillations')

	def PointOsc_SKTable(self):
		self.weightOsc_SKpaper = np.array([])
		units = nsq.Const()

		for k,(nu,E,cz,mod) in enumerate(zip(self.Ipnu, self.Enu, self.Cz, self.Mode)):
		# Get P_{x->ipnu} probabilities
			weight = 0.0
			if nu>0:
				rho = 0
			elif nu<0:
				rho = 1
			else:
				print('What?! No identified neutrino flavour')

			if cz>0:
				coszen = np.array([0.99999*cz,1.00001*cz])
			else:
				coszen = np.array([1.00001*cz,0.99999*cz])

			ener = np.array([0.99*E,1.01*E])
			AtmOsc = nsq.nuSQUIDSAtm(coszen,ener*units.GeV,3,nsq.NeutrinoType.both,False)
			AtmOsc.Set_rel_error(1.0e-5);
			AtmOsc.Set_abs_error(1.0e-5);
			AtmOsc.Set_MixingAngle(0,1,math.asin(math.sqrt(0.304)))
			AtmOsc.Set_MixingAngle(0,2,math.asin(math.sqrt(0.0219)))
			AtmOsc.Set_MixingAngle(1,2,math.asin(math.sqrt(0.5)))
			AtmOsc.Set_SquareMassDifference(1,7.53e-05)
			AtmOsc.Set_SquareMassDifference(2,0.0024)
			AtmOsc.Set_CPPhase(0,2,0)

			if abs(mod) < 30:
				for i in range(2):
					in_state = np.zeros((2,2,2,3))
					for ii in range(2):
						for ij in range(2):
							for ik in range(2):
								in_state[ii][ij][ik][i] = 1
					Ffactor = FluxFactor(i, nu, self.Flux_nue[k], self.Flux_nueb[k], self.Flux_numu[k], self.Flux_numub[k])
					AtmOsc.Set_initial_state(in_state,nsq.Basis.flavor)
					AtmOsc.EvolveState()
					j = int(abs(nu) / 2) % 6
					prob = 0
					for _ in range(20):
						prob += AtmOsc.EvalFlavor(j, cz, E*units.GeV, rho, True)
					prob /= 20
					weight += prob*Ffactor
			else:
				weight = 1.0

			self.weightOsc_SKpaper = np.append(self.weightOsc_SKpaper,weight)
			if math.isnan(weight): print(self.weightOsc_SKpaper[k],weight,cz,E)
		print('Done with oscillations')

	# def PointOsc_SKTable(self):
	# 	self.weightOsc_SKpaper = np.zeros_like(self.Enu)
	# 	units = nsq.Const()
	# 	interactions = False
		
	# 	AtmOsc = nsq.nuSQUIDSAtm(self.cth_nodes,self.energy_nodes*units.GeV,3,nsq.NeutrinoType.both,interactions)
	# 	AtmOsc.Set_rel_error(1.0e-4);
	# 	AtmOsc.Set_abs_error(1.0e-4);
	# 	AtmOsc.Set_MixingAngle(0,1,math.asin(math.sqrt(0.304)))
	# 	AtmOsc.Set_MixingAngle(0,2,math.asin(math.sqrt(0.0219)))
	# 	AtmOsc.Set_MixingAngle(1,2,math.asin(math.sqrt(0.5)))
	# 	AtmOsc.Set_SquareMassDifference(1,7.53e-05)
	# 	AtmOsc.Set_SquareMassDifference(2,0.0024)
	# 	AtmOsc.Set_CPPhase(0,2,0)
	# 	AtmOsc.Set_initial_state(self.InitialFlux,nsq.Basis.flavor)
	# 	AtmOsc.EvolveState()

	# 	neuflavor=0
	# 	for i,(E,cz) in enumerate(zip(self.Enu, self.Cz)):
	# 		if self.Ipnu[i] > 0 :
	# 			neutype = 0
	# 		else:
	# 			neutype = 1
	# 		if np.abs(self.Ipnu[i]) == 12:
	# 			neuflavor = 0
	# 		elif np.abs(self.Ipnu[i]) == 14:
	# 			neuflavor = 1
	# 		elif np.abs(self.Ipnu[i]) == 16:
	# 			neuflavor = 2
	# 		self.weightOsc_SKpaper[i] = AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype, True)
	# 	print('Done with oscillations')


	# def PointOsc_SKBest(self):
	# 	self.weightOsc_SKbest = np.zeros_like(self.Enu)
	# 	units = nsq.Const()
	# 	interactions = False
		
	# 	AtmOsc = nsq.nuSQUIDSAtm(self.cth_nodes,self.energy_nodes*units.GeV,3,nsq.NeutrinoType.both,interactions)
	# 	AtmOsc.Set_rel_error(1.0e-4);
	# 	AtmOsc.Set_abs_error(1.0e-4);
	# 	AtmOsc.Set_MixingAngle(0,1,math.asin(math.sqrt(0.304)))
	# 	AtmOsc.Set_MixingAngle(0,2,math.asin(math.sqrt(0.0219)))
	# 	AtmOsc.Set_MixingAngle(1,2,math.asin(math.sqrt(0.588)))
	# 	AtmOsc.Set_SquareMassDifference(1,7.53e-05)
	# 	AtmOsc.Set_SquareMassDifference(2,0.0025)
	# 	AtmOsc.Set_CPPhase(0,2,4.18)
	# 	AtmOsc.Set_initial_state(self.InitialFlux,nsq.Basis.flavor)
	# 	AtmOsc.EvolveState()

	# 	neuflavor=0
	# 	for i,(E,cz) in enumerate(zip(self.Enu, self.Cz)):
	# 		if self.Ipnu[i] > 0 :
	# 			neutype = 0
	# 		else:
	# 			neutype = 1
	# 		if np.abs(self.Ipnu[i]) == 12:
	# 			neuflavor = 0
	# 		elif np.abs(self.Ipnu[i]) == 14:
	# 			neuflavor = 1
	# 		elif np.abs(self.Ipnu[i]) == 16:
	# 			neuflavor = 2
	# 		self.weightOsc_SKbest[i] = AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype, True)
	# 		if i==1: 
	# 			print(AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype, True), self.Ipnu[i], neuflavor)
	# 			print(AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype), self.Ipnu[i], neuflavor)
	# 	print('Done with oscillations')


	# def InitialFlux(self):
	# 	flux = nuflux.makeFlux('IPhonda2014_sk_solmin')
	# 	E_min = 0.1
	# 	E_max = 4.0e2
	# 	E_nodes = 100
	# 	energy_range = nsq.logspace(E_min,E_max,E_nodes)
	# 	energy_nodes = nsq.logspace(E_min,E_max,E_nodes)
	# 	cth_min = -1.0
	# 	cth_max = 1.0
	# 	cth_nodes = 40
	# 	cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)
	# 	neutrino_flavors = 3

	# 	#Initialize the flux
	# 	AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
	# 	for ic,nu_cos_zenith in enumerate(cth_nodes):
	# 		for ie,nu_energy in enumerate(energy_range):
	# 			AtmInitialFlux[ic][ie][0][0] = 1 #flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
	# 			AtmInitialFlux[ic][ie][1][0] = 1 #flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
	# 			AtmInitialFlux[ic][ie][0][1] = 1 #flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
	# 			AtmInitialFlux[ic][ie][1][1] = 1 #flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
	# 			AtmInitialFlux[ic][ie][0][2] = 0.  # nutau
	# 			AtmInitialFlux[ic][ie][1][2] = 0.  # nutau bar
	# 	self.energy_nodes = energy_nodes
	# 	self.cth_nodes = cth_nodes
	# 	self.InitialFlux = AtmInitialFlux
