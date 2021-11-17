# Requirements
import os
import numpy as np
import argparse
import time
import h5py
# Local libraries
from TrueRing import TrueRingConstructor, numberOfNeutrons
from RecoRing import RecoRing
from RandomGenerator import RecoDists
from GenieHDF5 import FluxFactor, GenieSimulation
from pythiaDecay import pythiaDecay
from RecoNonFC import nonFCReco
from reWeight import simMatrix, SKMatrix, modeIndex

print('Super-Kamiokande atmospheric neutrino simulation')

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("in_hdf5filename", type=str, nargs='?', default='NULL')
	parser.add_argument("outfilename", type=str, nargs='?', default='NULL')
	parser.add_argument("-v", dest='verbo', default=False, action='store_true')
	args = parser.parse_args()
	verbose = args.verbo
	genie_input = args.in_hdf5filename
	output = args.outfilename
	if genie_input != 'NULL' and output == 'NULL':
		indir, fname = os.path.split(genie_input)
		output = 'data/output/sksim.'+fname

	# Read GENIE GST file in hdf5 format
	event = GenieSimulation(genie_input)

	# Read reco distributions
	rd = RecoDists()

	toposample = {'FC':0, 'PC-Stop':14, 'PC-Thru':15, 'UpMu-Stop':16, 'UpMu-Thru':17, 'UpMu-Show':18, 'None':-1}

	# Variable declaration for out file
	# True variables
	ipnu        = np.array([], dtype=np.double)
	pnu         = np.array([], dtype=np.double)
	dirnu_x     = np.array([], dtype=np.double)
	dirnu_y     = np.array([], dtype=np.double)
	dirnu_z     = np.array([], dtype=np.double)
	cz          = np.array([], dtype=np.double)
	azi         = np.array([], dtype=np.double)
	fluxho_numu = np.array([], dtype=np.double)
	fluxho_nue  = np.array([], dtype=np.double)
	fluxho_numub= np.array([], dtype=np.double)
	fluxho_nueb = np.array([], dtype=np.double)
	weightOsc_SKpaper   = np.array([], dtype=np.double)
	weightOsc_SKbest   = np.array([], dtype=np.double)
	weightSim   = np.array([], dtype=np.double)
	plep        = np.array([], dtype=np.double)
	dirlep_x    = np.array([], dtype=np.double)
	dirlep_y    = np.array([], dtype=np.double)
	dirlep_z    = np.array([], dtype=np.double)
	## Reco variables
	reco_pmax   = np.array([], dtype=np.double)
	evis        = np.array([], dtype=np.double)
	reco_dir_x  = np.array([], dtype=np.double)
	reco_dir_y  = np.array([], dtype=np.double)
	reco_dir_z  = np.array([], dtype=np.double)
	ip          = np.array([], dtype=np.double)
	nring       = np.array([], dtype=np.double)
	muedk       = np.array([], dtype=np.double)
	neutron     = np.array([], dtype=np.double)
	itype       = np.array([], dtype=np.double)
	imass       = np.array([], dtype=np.double)
	mode        = np.array([], dtype=np.double)

	# Enable decays of unstable particles
	pd = pythiaDecay()

	for i, nu in enumerate(event.Ipnu):
		if i>10000: break
		# if verbose:
		if i%1000 == 0:
			print('----------------------------------------------------')
			print('----------- Event number ',i,'----------------------')
			print('----------------------------------------------------')
	# Get variables from Genie Simulation
		if abs(event.Ipnu[i])==16 and event.NC[i]: # Removing tau NC from simulation
			continue
	# Neutrino
		Pnu_v = np.array([event.Pxnu[i], event.Pynu[i], event.Pznu[i]]) # Neutrino momentum vector
	# Primary lepton
		Plep_v = np.array([event.Pxlep[i], event.Pylep[i], event.Pzlep[i]]) # Prompt lepton momentum vector
	# Hadronic system arrangement
		P_fp_v = np.transpose(np.array([event.Pxhad[i], event.Pyhad[i], event.Pzhad[i]])) # Momenta
		P_fp = event.Phad[i] # Momentum
		pdg_fp = event.PDGhad[i] # PDG IDs
		E_fp = event.Ehad[i] # Energies
	# Merge all interaction products
		pdgs = np.append(event.PDGlep[i], pdg_fp)
		Es = np.append(event.Elep[i], E_fp)
		Ps = np.append(event.Plep[i], P_fp)
		Pvs = np.vstack((Plep_v, P_fp_v))
	# Anything to decay?
		decayed = np.array([], dtype='int32')
		for k, particle in enumerate(pdgs):
			if pd.canDecay(particle):
				decayed = np.append(decayed, k)
				daug_id, daug_E, daug_p, daug_pv = pd.decay(particle, Es[k], Pvs[k])
				if daug_id.size>0:
					pdgs = np.append(pdgs,daug_id)
					Es = np.append(Es,daug_E)
					Pvs = np.vstack((Pvs,daug_pv))
					Ps = np.append(Ps,daug_p)
		# Remove decayed particles
		pdgs = np.delete(pdgs, decayed)
		Es = np.delete(Es, decayed)
		Ps = np.delete(Ps, decayed)
		Pvs = np.delete(Pvs, decayed, 0)
	# True Ring constructor, baseline for reconstructing the event
		TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir = TrueRingConstructor(pdgs, Es, Ps, Pvs)
		if TrueNRing == 0:
			if verbose:
				print('Zero true Rings')
			continue

	# SK topology classification: Detailed reconstruction only for FCs, the rest are much simpler
		if event.TopologySample[i]=='FC':
		# Reconstructing rings
			RRing = RecoRing(TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir, rd, event.Mode[i])
			if RRing.NRing == 1:
				RRing.DecayE(event.Ipnu[i], event.CC[i], event.Mode[i])
			RRing.SKType(event.Ipnu[i], event.CC[i])
			if RRing.NRing < 1 or RRing.Type < 0:
				if verbose:
					print('Zero reconstructed Rings')
				del RRing
				continue
		# Verbosity
			if verbose:
				print('Neutrino direction: ', Pnu_v / event.Enu[i])
				print('Neutrino energy: ', event.Enu[i])
				print('Neutrino id: ', event.Ipnu[i])
				print('Event with ', RRing.NRing, ' reco rings')
				print('Visible energy is ', RRing.Evis, ' GeV')
				print('Event type is ', RRing.Type)
				print('Reconstructed total direction: ', RRing.TotDir)
				print('---------------------------------')

		elif toposample[event.TopologySample[i]]>1:
			nonfcType = toposample[event.TopologySample[i]]
			if TrueNRing == 1:
				RnonFC = nonFCReco(TrueRingPDG[0], nonfcType, rd, TrueRingDir, TrueRingP[0])
			else:
				maxP = np.argmax(TrueRingP)
				RnonFC = nonFCReco(TrueRingPDG[maxP], nonfcType, rd, TrueRingDir[maxP], TrueRingP[maxP])
			nonfcEvis = RnonFC.Momentum
			nonfcTotXDir = RnonFC.Direction[0]
			nonfcTotYDir = RnonFC.Direction[1]
			nonfcTotZDir = RnonFC.Direction[2]
		# Verbosity
			if verbose:
				print('Neutrino direction: ', Pnu_v / event.Enu[i])
				print('Neutrino energy: ', event.Enu[i])
				print('Neutrino id: ', event.Ipnu[i])
				print('Visible energy is ', RnonFC.Momentum, ' GeV')
				print('Event type is ', nonfcType)
				print('Reconstructed total direction: ', RnonFC.Direction)
				print('---------------------------------')
		else:
			continue

	# Fill true variables
		ipnu       = np.append(ipnu, event.Ipnu[i])
		pnu        = np.append(pnu, event.Enu[i])
		dirnu_x    = np.append(dirnu_x, Pnu_v[0] / event.Enu[i])
		dirnu_y    = np.append(dirnu_y, Pnu_v[1] / event.Enu[i])
		dirnu_z    = np.append(dirnu_z, Pnu_v[2] / event.Enu[i])
		azi        = np.append(azi, event.Azi[i])
		cz         = np.append(cz, event.Cz[i])
		fluxho_nue  = np.append(fluxho_nue, event.Flux_nue[i])
		fluxho_nueb  = np.append(fluxho_nueb, event.Flux_nueb[i])
		fluxho_numu  = np.append(fluxho_numu, event.Flux_numu[i])
		fluxho_numub  = np.append(fluxho_numub, event.Flux_numub[i])
		weightSim  = np.append(weightSim, event.FluxWeight[i])
		weightOsc_SKpaper = np.append(weightOsc_SKpaper, event.weightOsc_SKpaper[i])
		weightOsc_SKbest = np.append(weightOsc_SKbest, event.weightOsc_SKbest[i])
		plep       = np.append(plep, event.Plep[i])
		dirlep_x   = np.append(dirlep_x, event.Pxlep[i] / event.Plep[i])
		dirlep_y   = np.append(dirlep_y, event.Pylep[i] / event.Plep[i])
		dirlep_z   = np.append(dirlep_z, event.Pzlep[i] / event.Plep[i])
		mode       = np.append(mode, event.Mode[i])
	# Fill reco variables
		if event.TopologySample[i]=='FC':
			reco_pmax  = np.append(reco_pmax, RRing.MERMomentum)
			evis       = np.append(evis, RRing.Evis)
			reco_dir_x = np.append(reco_dir_x, RRing.TotDir[0])
			reco_dir_y = np.append(reco_dir_y, RRing.TotDir[1])
			reco_dir_z = np.append(reco_dir_z, RRing.TotDir[2])
			ip         = np.append(ip, RRing.MERIP)
			nring      = np.append(nring, RRing.NRing)
			muedk      = np.append(muedk, RRing.MuEdk)
			neutron    = np.append(neutron, numberOfNeutrons(pdgs))
			itype      = np.append(itype, RRing.Type)
			imass      = np.append(imass, RRing.Imass)
			del RRing
		else:
			reco_pmax  = np.append(reco_pmax, -9999.)
			evis       = np.append(evis, nonfcEvis)
			reco_dir_x = np.append(reco_dir_x, nonfcTotXDir)
			reco_dir_y = np.append(reco_dir_y, nonfcTotYDir)
			reco_dir_z = np.append(reco_dir_z, nonfcTotZDir)
			ip         = np.append(ip, 0)
			nring      = np.append(nring, 0)
			muedk      = np.append(muedk, 0)
			neutron    = np.append(neutron, 0)
			itype      = np.append(itype, nonfcType)
			imass      = np.append(imass, -9999.)

	# Applying weights to match SK's public event rate tables
	W = simMatrix(itype, ipnu, mode, weightOsc_SKpaper) # Rate matrix from this simulation
	W0= SKMatrix() # Rate matrix from SK's paper
	weightReco = np.zeros(np.size(weightOsc_SKpaper))
	for i,t in enumerate(itype):
		ti = int(t)
		j = modeIndex(ipnu[i], mode[i])
		if ti>-1 and ti<16:
			# print(ipnu[i], mode[i], j)
			if W[ti][j]>0:
				weightReco[i] = W0[ti][j] / W[ti][j]
			else:
				print('Simulation weight is zero!')
				weightReco[i] = 0
		else:
			weightReco[i] = 0

	# Saving data
	with h5py.File(output, 'w') as hf:
		hf.create_dataset('ipnu', data=ipnu, compression='gzip')
		hf.create_dataset('pnu', data=pnu, compression='gzip')
		hf.create_dataset('dirnuX', data=dirnu_x, compression='gzip')
		hf.create_dataset('dirnuY', data=dirnu_y, compression='gzip')
		hf.create_dataset('dirnuZ', data=dirnu_z, compression='gzip')
		hf.create_dataset('azi', data=azi, compression='gzip')
		hf.create_dataset('cz', data=cz, compression='gzip')
		hf.create_dataset('fluxho_nue', data=fluxho_nue, compression='gzip')
		hf.create_dataset('fluxho_nueb', data=fluxho_nueb, compression='gzip')
		hf.create_dataset('fluxho_numu', data=fluxho_numu, compression='gzip')
		hf.create_dataset('fluxho_numub', data=fluxho_numub, compression='gzip')
		hf.create_dataset('plep', data=plep, compression='gzip')
		hf.create_dataset('dirlepX', data=dirlep_x, compression='gzip')
		hf.create_dataset('dirlepY', data=dirlep_y, compression='gzip')
		hf.create_dataset('dirlepZ', data=dirlep_z, compression='gzip')
		hf.create_dataset('pmax', data=reco_pmax, compression='gzip')
		hf.create_dataset('evis', data=evis, compression='gzip')
		hf.create_dataset('recodirX', data=reco_dir_x, compression='gzip')
		hf.create_dataset('recodirY', data=reco_dir_y, compression='gzip')
		hf.create_dataset('recodirZ', data=reco_dir_z, compression='gzip')
		hf.create_dataset('ip', data=ip, compression='gzip')
		hf.create_dataset('nring', data=nring, compression='gzip')
		hf.create_dataset('muedk', data=muedk, compression='gzip')
		hf.create_dataset('neutron', data=neutron, compression='gzip')
		hf.create_dataset('itype', data=itype, compression='gzip')
		hf.create_dataset('imass', data=imass, compression='gzip')
		hf.create_dataset('mode', data=mode, compression='gzip')
		hf.create_dataset('weightSim', data=weightSim, compression='gzip')
		hf.create_dataset('weightOsc_SKpaper', data=weightOsc_SKpaper, compression='gzip')
		hf.create_dataset('weightOsc_SKbest', data=weightOsc_SKbest, compression='gzip')
		hf.create_dataset('weightReco', data=weightReco, compression='gzip')

if __name__ == '__main__':
	main()
