# Requirements
import numpy as np
import argparse
import time
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
	parser.add_argument("in_hdf5filename", type=str, nargs='?', default='/home/pablofer/old.AtmNuCombination/utils/atm_genie.root.hdf5')
	# parser.add_argument("in_hdf5filename", type=str, nargs='?', default='/home/pablofer/IceCube-SuperK_AtmNu/SuperK/AtmNu_MC/atm_tau/gntp.0.gst.root.hdf5')
	# parser.add_argument("in_hdf5filename", type=str, nargs='?', default='/home/pablofer/IceCube-SuperK_AtmNu/SuperK/AtmNu_MC/genie_atm_full_1p5M.root.hdf5')
	parser.add_argument("outfilename", type=str, nargs='?', default='data/testfcmc.hdf5')
	parser.add_argument("-v", dest='verbo', default=False, action='store_true')
	args = parser.parse_args()
	verbose = args.verbo
	genie_input = args.in_hdf5filename
	output = args.outfilename

	# Read GENIE GST file in hdf5 format
	event = GenieSimulation(genie_input)

	# Read reco distributions
	# start_time = time.time()
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
	weightOsc   = np.array([], dtype=np.double)
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

	# Start loop over Genie simulated
	start_time = time.time()

	for i, nu in enumerate(event.Ipnu):
		# if i>1000: break
		# if verbose:
		if i%1000 == 0:
			print('----------------------------------------------------')
			print('----------- Event number ',i,'----------------------')
			print('----------------------------------------------------')

	# Get variables from Genie Simulation
		if abs(event.Ipnu[i])==16 and event.NC[i]:
			continue
		# 	# print('Removing tau NC...')
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
		print('pdgs before:', pdgs)
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
		print('pdgs after:', pdgs)
		Es = np.delete(Es, decayed)
		Ps = np.delete(Ps, decayed)
		Pvs = np.delete(Pvs, decayed, 0)
			
	# True Ring constructor, baseline for reconstructing the event
		TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir = TrueRingConstructor(pdgs, Es, Ps, Pvs)
		if TrueNRing == 0:
			if verbose:
				print('Zero Rings')
			continue

	# SK topology classification: Detailed reconstruction only for FCs, the rest are much simpler
		if event.TopologySample[i]=='FC':

		# Reconstructing rings
			RRing = RecoRing(TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir, rd, event.Mode[i])
			if RRing.NRing == 1:
				RRing.DecayE(event.Ipnu[i], event.CC[i], event.Mode[i])
			RRing.SKType(event.Ipnu[i], event.CC[i])
			if RRing.NRing < 1 or RRing.Type < 0:
				del RRing
				continue
			# # Fine tuning
			# RRing.mendSGE()
			# RRing.mendSGM(event.Ipnu[i])
			# RRing.mendMG(event.Ipnu[i])
			# RRing.mendMR(event.Ipnu[i])

		# Verbosity
			if verbose:
				if RRing.NRing == 1:
					print('Single Ring event')
					print('Visible energy is ', RRing.Evis, ' GeV')
					print('Event IP is ', RRing.IP)
					print('Number of Michel electrons: ', RRing.MuEdk)
				if RRing.NRing > 1:
					print('Multi Ring event')
					print('Event with ', RRing.NRing, ' reco rings')
					print('Visible energy is ', RRing.Evis, ' GeV')
					print('Event IP is ', RRing.IP)
					print('Event MER IP is ', RRing.MERIP)
					print('Event type is ', RRing.Type)

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
		
		else:
			continue

	# Verbosity
		if verbose:
			print('Event particle IDs', pdgs)
			print('Particle energies', Es)
			print('Particle momenta', Ps)
			print('Neutrino energy: ',event.Enu[i],' GeV')
			print('Neutrino flavour: ',event.Ipnu[i])
			print('Lepton ID: ',event.PDGlep)
			print('Lepton momentum: ',event.LeptonMomentum, 'GeV')
			if abs(event.Mode[i])<30:
				print('Interacting CC with code: ', event.Mode[i])
			else:
				print('Interacting NC with code: ', event.Mode[i])
			print('To be reconstructed as ',event.TopologySample)


	# Fill variables
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
		weightOsc = np.append(weightOsc, event.weightOsc[i])
		plep       = np.append(plep, event.Plep[i])
		dirlep_x   = np.append(dirlep_x, event.Pxlep[i] / event.Plep[i])
		dirlep_y   = np.append(dirlep_y, event.Pylep[i] / event.Plep[i])
		dirlep_z   = np.append(dirlep_z, event.Pzlep[i] / event.Plep[i])
		mode       = np.append(mode, event.Mode[i])
		# reco variables
		if event.TopologySample[i]=='FC':
			# print(event.weightOsc[i], Pnu_v[2] / event.Enu[i], RRing.TotDir[2])
			# print(event.Pzlep[i] / event.Plep[i])
			# if np.dot(RRing.TotDir,Pnu_v / event.Enu[i])<0.5:
			# 	print('TrueRingPDG:', TrueRingPDG)
			# 	print('dirnu:', Pnu_v / event.Enu[i])
			# 	print('RRing.TotDir', RRing.TotDir)
			# 	print('TrueRingDir:',TrueRingDir)
			# 	print('RRing.Direction', RRing.Direction)
			# 	print('================================')
			print('================================')

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
	# W = simMatrix(itype, ipnu, mode, weightOsc*weightSim) # Rate matrix from this simulation
	W = simMatrix(itype, ipnu, mode, weightOsc) # Rate matrix from this simulation
	W0= SKMatrix() # Rate matrix from SK's paper
	weightReco = np.zeros(np.size(weightOsc))
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
		# print(weightOsc[i], weightReco[i], dirnu_z[i], reco_dir_z[i])
	# ww = np.multiply(weightOsc,weightReco)

	# Save data
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
		hf.create_dataset('weightOsc', data=weightOsc, compression='gzip')
		hf.create_dataset('itype', data=itype, compression='gzip')
		hf.create_dataset('imass', data=imass, compression='gzip')
		hf.create_dataset('mode', data=mode, compression='gzip')
		hf.create_dataset('weightReco', data=weightReco, compression='gzip')

if __name__ == '__main__':
	main()