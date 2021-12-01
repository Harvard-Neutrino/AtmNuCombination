import numpy as np
import pandas as pd
import h5py
import nuflux
import math
import nuSQuIDS as nsq
import nuSQUIDSTools
from math import asin, sqrt

class Reader:
	def __init__(self, experiment, filename):

		self.Experiment = experiment
		# if filename[-4:-1] == 'hdf5' or filename[-2:-1] == 'h5':
		if self.Experiment == 'Super-Kamiokande' or self.Experiment == 'SK' or self.Experiment == 'Super-Kamiokande':
			print(f'Processing simulation of {self.Experiment} experiment.')
			with h5py.File(filename,'r') as hf:
				d_evis = np.array(hf['evis'][()])
				d_recocz = np.array(hf['recodirZ'][()])
				d_truecz = np.array(hf['dirnuZ'][()])
				d_azi = np.array(hf['azi'][()])
				d_mode = np.array(hf['mode'][()], dtype=int)
				d_ipnu = np.array(hf['ipnu'][()], dtype=int)
				d_pnu = np.array(hf['pnu'][()])
				d_weightReco = np.array(hf['weightReco'][()])
				d_weightSim = np.array(hf['weightSim'][()])
				d_itype = np.array(hf['itype'][()], dtype=int)
			# To be removed at some point
			condition1 = (d_itype<14) * (d_itype>-1) 
			condition2 = (d_itype<16) * (d_itype>13) * (d_evis>1)
			condition = (condition2 + condition1) * (d_evis<400) * (d_pnu>0.1) * (d_pnu<400)
			self.EReco = d_evis[condition]
			self.CosZReco = d_recocz[condition]
			self.CosZTrue = d_truecz[condition]
			self.AziTrue = d_azi[condition]
			self.CC = np.abs(d_mode[condition]) < 30
			self.nuPDG = d_ipnu[condition]
			self.ETrue = d_pnu[condition]
			self.Weight = d_weightReco[condition] * d_weightSim[condition]
			self.Sample = d_itype[condition]
			self.Norm = 0.05411673990568942
			self.NumberOfSamples = 16
			self.Erec_min = 0.1
			self.NumberOfEvents = self.nuPDG.size
		elif self.Experiment == 'IceCube' or self.Experiment == 'IC' or self.Experiment == 'DeepCore':
			print(f'Processing simulation of {self.Experiment} experiment.')
			input_data = pd.read_csv(filename)
			Time = 3.0*365*24*60*60
			meter_to_cm_sq = 1e4
			self.Erec_min = 1
			d_EReco = input_data['reco_energy'].to_numpy()
			d_CosZReco = np.cos(input_data['reco_zenith'].to_numpy())
			d_CosZTrue = np.cos(input_data['true_zenith'].to_numpy())
			d_AziTrue = np.cos(input_data['true_azimuth'].to_numpy())
			d_CC = input_data['current_type'].to_numpy()
			d_nuPDG = np.int_(input_data['pdg'].to_numpy())
			d_ETrue = input_data['true_energy'].to_numpy()
			d_Weight = input_data['weight'].to_numpy()
			d_Sample = np.int_(input_data['pid'].to_numpy())
			condition = (d_ETrue > 1) * (d_ETrue < 1e3) * (d_EReco > 1)
			self.EReco = d_EReco[condition]
			self.CosZReco = d_CosZReco[condition]
			self.CosZTrue = d_CosZTrue[condition]
			self.AziTrue = d_AziTrue[condition]
			self.CC = d_CC[condition]
			self.nuPDG = d_nuPDG[condition]
			self.ETrue = d_ETrue[condition]
			self.Weight = d_Weight[condition]
			self.Sample = d_Sample[condition]
			self.Norm = Time * meter_to_cm_sq
			self.NumberOfSamples = 2
			self.NumberOfEvents = self.nuPDG.size
		# print(self.NumberOfEvents)

	def SetOutFile(self,o):
		self.outfile = o

	def Binning(self):
		if self.Experiment == 'Super-Kamiokande' or self.Experiment == 'SK':
			sge_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
			sgm_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
			sgsrpi0ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 500])
			sgmrpi0ebins = np.array([0.1, 0.15, 0.25, 0.4, 0.63, 500])
			mge_ebins = np.array([1.3, 2.5, 5., 10., 500.])
			mgm_ebins = np.array([1.3, 3.0, 500.])
			mre_ebins = np.array([1.3, 2.5, 5.0, 500.])
			mrm_ebins = np.array([0.6, 1.3, 2.5, 5., 500.])
			mro_ebins = np.array([1.3, 2.5, 5.0, 10., 500.])
			pcs_ebins = np.array([0.1, 20., 1.0e5])
			pct_ebins = np.array([0.1, 10.0, 50., 1.0e3, 1.0e5])
			z10bins = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
			z1bins = np.array([-1, 1.0])
			self.EnergyBins = {0:sge_ebins, 1:sge_ebins, 2:sgsrpi0ebins, 3:sgm_ebins, 4:sgm_ebins, 5:sgm_ebins, 6:sgmrpi0ebins,
			7:mge_ebins, 8:mge_ebins, 9:mgm_ebins, 10:mre_ebins, 11:mre_ebins, 12:mrm_ebins, 13:mro_ebins, 14:pcs_ebins, 15:pct_ebins}
			self.CzBins = {0:z10bins, 1:z1bins, 2:z1bins, 3:z10bins, 4:z10bins, 5:z1bins, 6:z1bins,
			7:z10bins, 8:z10bins, 9:z10bins, 10:z10bins, 11:z10bins, 12:z10bins, 13:z10bins, 14:z10bins, 15:z10bins}
			self.MaxNumberOfEnergyBins = 5
			self.MaxNumberOfCzBins = 10
		elif self.Experiment == 'IceCube' or self.Experiment == 'IC' or self.Experiment == 'DeepCore':
			Erec_max = 1e4
			NErec = 40
			erec = np.logspace(np.log10(self.Erec_min), np.log10(Erec_max), NErec+1, endpoint = True)
			z10bins = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
			self.EnergyBins = {0:erec, 1:erec}
			self.CzBins = {0:z10bins, 1:z10bins}
			self.MaxNumberOfEnergyBins = 40
			self.MaxNumberOfCzBins = 10


	def BFOscillator(self,th23,dm31,dcp=4.01425,th13=0.14957):#  230*Radian 8.57 * Radian
		units = nsq.Const()
		Radian = np.pi / 180.
		interactions = False
		neutrino_flavors = 3
		AtmOsc = nsq.nuSQUIDSAtm(self.cth_nodes,self.energy_nodes*units.GeV,neutrino_flavors,nsq.NeutrinoType.both,interactions)
		AtmOsc.Set_rel_error(1.0e-4);
		AtmOsc.Set_abs_error(1.0e-4);
		AtmOsc.Set_MixingAngle(0,1, 33.44 * Radian)
		AtmOsc.Set_MixingAngle(0,2, th13)
		AtmOsc.Set_SquareMassDifference(1,7.42e-5)
		AtmOsc.Set_MixingAngle(1,2, th23)
		AtmOsc.Set_SquareMassDifference(2,dm31)
		AtmOsc.Set_CPPhase(0,2,dcp)
		AtmOsc.Set_initial_state(self.AtmInitialFlux,nsq.Basis.flavor)
		AtmOsc.EvolveState()

		self.weightOscBF = np.zeros(self.NumberOfEvents)
		neuflavor=0
		for i,(E,cz) in enumerate(zip(self.ETrue, self.CosZTrue)):
			if self.nuPDG[i] > 0 :
				neutype = 0
			else:
				neutype = 1
			if np.abs(self.nuPDG[i]) == 12:
				neuflavor = 0
			elif np.abs(self.nuPDG[i]) == 14:
				neuflavor = 1
			elif np.abs(self.nuPDG[i]) == 16:
				neuflavor = 2
			self.weightOscBF[i] = AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype)

	def BFOscillatorxml(self,neutrino_flavors, Sin2Theta12=0, Sin2Theta13=0, Sin2Theta23=0, Dm221=0, Dm231=0, dCP=0, Ordering='normal'):
		units = nsq.Const()
		interactions = False
		AtmOsc = nsq.nuSQUIDSAtm(self.cth_nodes,self.energy_nodes*units.GeV,neutrino_flavors,nsq.NeutrinoType.both,interactions)
		AtmOsc.Set_rel_error(1.0e-4);
		AtmOsc.Set_abs_error(1.0e-4);
		AtmOsc.Set_MixingAngle(0,1, asin(sqrt(Sin2Theta12)))
		AtmOsc.Set_MixingAngle(0,2, asin(sqrt(Sin2Theta13)))
		AtmOsc.Set_SquareMassDifference(1, Dm221)
		AtmOsc.Set_MixingAngle(1,2, asin(sqrt(Sin2Theta23)))
		if Ordering=='normal':
			AtmOsc.Set_SquareMassDifference(2,Dm231)
		else:
			AtmOsc.Set_SquareMassDifference(2,Dm221-Dm231)
		AtmOsc.Set_CPPhase(0,2,dCP)
		AtmOsc.Set_initial_state(self.AtmInitialFlux,nsq.Basis.flavor)
		AtmOsc.EvolveState()

		self.weightOscBF = np.zeros(self.NumberOfEvents)
		neuflavor=0
		for i,(E,cz) in enumerate(zip(self.ETrue, self.CosZTrue)):
			if self.nuPDG[i] > 0 :
				neutype = 0
			else:
				neutype = 1
			if np.abs(self.nuPDG[i]) == 12:
				neuflavor = 0
			elif np.abs(self.nuPDG[i]) == 14:
				neuflavor = 1
			elif np.abs(self.nuPDG[i]) == 16:
				neuflavor = 2
			self.weightOscBF[i] = AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype)

	def GridOscPar(self,nt23,nm31,ncp):
		return np.zeros((nt23,nm31,ncp,self.NumberOfEvents))
		# self.weightOsc = np.zeros((nt23+1,nm31+1,ncp+1,self.NumberOfEvents))

	def OscillatorM(self, j, th23, dm31, dcp, th13, IO, return_dict):
		units = nsq.Const()
		Radian = np.pi / 180.
		interactions = False
		neutrino_flavors = 3
		AtmOsc = nsq.nuSQUIDSAtm(self.cth_nodes,self.energy_nodes*units.GeV,neutrino_flavors,nsq.NeutrinoType.both,interactions)
		AtmOsc.Set_rel_error(1.0e-4);
		AtmOsc.Set_abs_error(1.0e-4);
		AtmOsc.Set_MixingAngle(0,1, 33.44 * Radian)
		AtmOsc.Set_MixingAngle(0,2, th13)
		AtmOsc.Set_SquareMassDifference(1,7.42e-5)
		AtmOsc.Set_MixingAngle(1,2, th23)
		if IO==1:
			AtmOsc.Set_SquareMassDifference(2,7.42e-5-dm31)
		else:
			AtmOsc.Set_SquareMassDifference(2,dm31)
		AtmOsc.Set_CPPhase(0,2,dcp)
		AtmOsc.Set_initial_state(self.AtmInitialFlux,nsq.Basis.flavor)
		AtmOsc.EvolveState()

		w = np.zeros(self.NumberOfEvents)

		for i,(E,cz) in enumerate(zip(self.ETrue, self.CosZTrue)):
			if self.nuPDG[i] > 0 :
				neutype = 0
			else:
				neutype = 1
			if np.abs(self.nuPDG[i]) == 12:
				neuflavor = 0
			elif np.abs(self.nuPDG[i]) == 14:
				neuflavor = 1
			elif np.abs(self.nuPDG[i]) == 16:
				neuflavor = 2
			w[i] = AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype)

		return_dict[j] = w

	def OscillatorMxml(self, neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering='normal'):
		units = nsq.Const()
		interactions = False
		AtmOsc = nsq.nuSQUIDSAtm(self.cth_nodes,self.energy_nodes*units.GeV,neutrino_flavors,nsq.NeutrinoType.both,interactions)
		AtmOsc.Set_rel_error(1.0e-4);
		AtmOsc.Set_abs_error(1.0e-4);
		AtmOsc.Set_MixingAngle(0,1, asin(sqrt(t12)))
		AtmOsc.Set_MixingAngle(0,2, asin(sqrt(t13)))
		AtmOsc.Set_SquareMassDifference(1,dm21)
		AtmOsc.Set_MixingAngle(1,2, asin(sqrt(t23)))
		if Ordering=='normal':
			AtmOsc.Set_SquareMassDifference(2,dm31)
		else:
			AtmOsc.Set_SquareMassDifference(2,dm21-dm31)
		AtmOsc.Set_CPPhase(0,2,dcp)
		AtmOsc.Set_initial_state(self.AtmInitialFlux,nsq.Basis.flavor)
		AtmOsc.EvolveState()

		w = np.zeros(self.NumberOfEvents)

		for i,(E,cz) in enumerate(zip(self.ETrue, self.CosZTrue)):
			if self.nuPDG[i] > 0 :
				neutype = 0
			else:
				neutype = 1
			if np.abs(self.nuPDG[i]) == 12:
				neuflavor = 0
			elif np.abs(self.nuPDG[i]) == 14:
				neuflavor = 1
			elif np.abs(self.nuPDG[i]) == 16:
				neuflavor = 2
			w[i] = AtmOsc.EvalFlavor(neuflavor, cz, E*units.GeV, neutype)

		return w

	def Chi2Calculator(self, neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering):
		wOsc = self.OscillatorMxml(neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering)
		X2=0
		for s in range(self.NumberOfSamples):
			wBF = self.weightOscBF[self.Sample==s] * self.Weight[self.Sample==s]
			Exp, dx, dy = np.histogram2d(self.EReco[self.Sample==s], self.CosZReco[self.Sample==s], bins=(self.EnergyBins[s], self.CzBins[s]), weights=wBF*self.Norm)
			w = wOsc[self.Sample==s] * self.Weight[self.Sample==s]
			Obs, dx, dy = np.histogram2d(self.EReco[self.Sample==s], self.CosZReco[self.Sample==s], bins=(self.EnergyBins[s], self.CzBins[s]), weights=w*self.Norm)
			for O,E in zip(np.ravel(Obs),np.ravel(Exp)):
				if O==0 or E==0: 
					pass
				else:
					X2 = X2 + 2 * (E - O + O * math.log(O/E))

		# print(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {X2}\n')
		with open(self.outfile,'a') as f:
			f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {X2}\n')
			f.flush()

		return X2

	def DefinSystematics(self):
		self.systBF = np.array([1.,1.,1.,1.,1.,1.])
		self.systSig = np.array([1.,1.,1.,1.,1.,1.])

	def Chi2CalculatorWsyst(self, syst, neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering):
		wOsc = self.OscillatorMxml(neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering)
		X2=0

		# norm = syst[0]
		# nnbar = syst[1]
		# eovermu = syst[2]
		# tilt = syst[3]
		# up = syst[4]
		# down = syst[5]

		systWeight = np.ones(self.NumberOfEvents) 
		tilt = (self.ETrue / self.E0Gam)**syst[3]
		norm = syst[0]
		nnbar = np.ones(self.NumberOfEvents)
		nnbar = nnbar[self.nuPDG<0] * syst[1]
		eovermu = np.ones(self.NumberOfEvents)
		eovermu = eovermu[abs(self.nuPDG)==12] * syst[2]
		zenith = np.ones(self.NumberOfEvents) 
		zenith = zenith[self.CosZTrue<0] - syst[4] * np.tanh(self.CosZTrue[self.CosZTrue<0])**2
		zenith = zenith[self.CosZTrue>=0] - syst[5] * np.tanh(self.CosZTrue[self.CosZTrue>=0])**2

		systWeight = systWeight * tilt * norm * nnbar * eovermu * zenith

		for s in range(self.NumberOfSamples):
			wBF = self.weightOscBF[self.Sample==s] * self.Weight[self.Sample==s]
			Exp, dx, dy = np.histogram2d(self.EReco[self.Sample==s], self.CosZReco[self.Sample==s], bins=(self.EnergyBins[s], self.CzBins[s]), weights=wBF*self.Norm)
			w = wOsc[self.Sample==s] * self.Weight[self.Sample==s] * systWeight[self.Sample==s]
			Obs, dx, dy = np.histogram2d(self.EReco[self.Sample==s], self.CosZReco[self.Sample==s], bins=(self.EnergyBins[s], self.CzBins[s]), weights=w*self.Norm)
			for O,E in zip(np.ravel(Obs),np.ravel(Exp)):
				if O==0 or E==0: 
					pass
				else:
					X2 = X2 + 2 * (E - O + O * math.log(O/E))

		for i,s in enumerate(syst):
			X2 = X2 + ((s-self.systBF[i]) / self.systSig[i])**2

		# print(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {X2}\n')
		# with open(self.outfile,'a') as f:
		# 	f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {X2}\n')
		# 	f.flush()

		return X2

		# minimize(fun, x0, args=(a,),

	def InitialFlux(self):
		if self.Experiment == 'Super-Kamiokande':
			flux = nuflux.makeFlux('IPhonda2014_sk_solmin')
			E_min = 0.1
			E_max = 4.0e2
		else:
			flux = nuflux.makeFlux('IPhonda2014_spl_solmin')
			E_min = 1.0
			E_max = 1.0e3

		E_nodes = 100
		energy_range = nsq.logspace(E_min,E_max,E_nodes)
		energy_nodes = nsq.logspace(E_min,E_max,E_nodes)
		cth_min = -1.0
		cth_max = 1.0
		cth_nodes = 40
		cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)
		neutrino_flavors = 3

		#Initialize the flux
		AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
		AtmIninue = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
		AtmIninueBar = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
		AtmIninum = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
		AtmIninumBar = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
		for ic,nu_cos_zenith in enumerate(cth_nodes):
			for ie,nu_energy in enumerate(energy_range):
				AtmInitialFlux[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
				AtmInitialFlux[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
				AtmInitialFlux[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
				AtmInitialFlux[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
				AtmInitialFlux[ic][ie][0][2] = 0.  # nutau
				AtmInitialFlux[ic][ie][1][2] = 0.  # nutau bar
		self.energy_nodes = energy_nodes
		self.cth_nodes = cth_nodes
		self.AtmInitialFlux = AtmInitialFlux
		
def OscPar(ndm31, dm31_min, dm31_max, nth23, th23_min, th23_max, ndcp, dcp_min, dcp_max):
	SqT23 = np.linspace(th23_min, th23_max, nth23, endpoint = True)
	T23 = np.arcsin(np.sqrt(SqT23))
	DM31 = np.linspace(dm31_min, dm31_max, ndm31, endpoint = True)
	DCP = np.linspace(dcp_min, dcp_max, ndcp, endpoint = True)
	return DM31, T23, SqT23, DCP

		


''' Probably outdated

	def EventMatrix(self, ndm31, nt23):
		self.Events = np.zeros(self.NumberOfSamples, ndm31, nt23, self.MaxNumberOfEnergyBins, self.MaxNumberOfCzBins)
		self.EventsBestFit = np.zeros(self.NumberOfSamples, self.MaxNumberOfEnergyBins, self.MaxNumberOfCzBins)
		self.WE = np.zeros((self.NumberOfSamples, self.NumberOfEvents))
		self.WEBAR = np.zeros((self.NumberOfSamples, self.NumberOfEvents))
		self.WM = np.zeros((self.NumberOfSamples, self.NumberOfEvents))
		self.WMBAR = np.zeros((self.NumberOfSamples, self.NumberOfEvents))

	def EventBin(self):
		# Asign a bin number to each event
		dBin = np.zeros((2,self.NumberOfEvents), dtype=int)
		for i in range(self.NumberOfEvents):
			h, bins = np.histogram(self.EReco[i],self.EnergyBins[self.Sample[i]])
			index = np.where(h==1)
			dBin[0][i] = int(index[0][0])
			h, bins = np.histogram(self.CosZReco[i],self.CzBins[self.Sample[i]])
			index = np.where(h==1)
			dBin[1][i] = int(index[0][0])
		return dBin
'''




