# Requirements
import numpy as np
import argparse
# Local dependencies
from TrueRing import *
from RecoRing import *
from RandomGenerator import RecoDists
from GenieHDF5 import *
from pythiaDecay import *

import time

parser = argparse.ArgumentParser()
parser.add_argument("in_hdf5filename", type=str, nargs='?', default='../../../utils/atm_genie.root.hdf5')
parser.add_argument("outfilename", type=str, nargs='?', default='data/testfcmc.hdf5')
parser.add_argument("-v", dest='verbo', default=False, action='store_true')

args = parser.parse_args()
verbose = args.verbo
genie_input = args.in_hdf5filename
output = args.outfilename

# Read GENIE GST file
genie_input = '/home/pablofer/old.AtmNuCombination/utils/atm_genie.root.hdf5'
start_time = time.time()
event = GenieSimulation(genie_input)
print('Event loading:', time.time() - start_time, 'seconds')

# Read reco distributions
start_time = time.time()
rd = RecoDists()
print('Reco. distros.:', time.time() - start_time, 'seconds')

# Variable declaration for out tree
# True variables
ipnu       = np.array([], dtype=np.double)
pnu        = np.array([], dtype=np.double)
dirnu_x    = np.array([], dtype=np.double)
dirnu_y    = np.array([], dtype=np.double)
dirnu_z    = np.array([], dtype=np.double)
cz         = np.array([], dtype=np.double)
azi        = np.array([], dtype=np.double)
fluxho_numu = np.array([], dtype=np.double)
fluxho_nue  = np.array([], dtype=np.double)
fluxho_numub= np.array([], dtype=np.double)
fluxho_nueb = np.array([], dtype=np.double)
oscw       = np.array([], dtype=np.double)
plep       = np.array([], dtype=np.double)
dirlep_x   = np.array([], dtype=np.double)
dirlep_y   = np.array([], dtype=np.double)
dirlep_z   = np.array([], dtype=np.double)
## Reco variables
reco_pmax  = np.array([], dtype=np.double)
evis       = np.array([], dtype=np.double)
reco_dir_x = np.array([], dtype=np.double)
reco_dir_y = np.array([], dtype=np.double)
reco_dir_z = np.array([], dtype=np.double)
ip         = np.array([], dtype=np.double)
nring      = np.array([], dtype=np.double)
muedk      = np.array([], dtype=np.double)
neutron    = np.array([], dtype=np.double)
itype      = np.array([], dtype=np.double)
mode       = np.array([], dtype=np.double)

# Enable decays of unstable particles
pd = pythiaDecay()

# Start loop over Genie simulated
start_time = time.time()

for i, nu in enumerate(event.Ipnu):

	# if i>10000: continue
	
	# if verbose:
	print('----------------------------------------------------')
	print('----------- Event number ',i,'--------------------')
	print('----------------------------------------------------')

# Get variables from Genie Simulation
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

	# print(pdgs)
	# print(Es)
	# print(Ps)

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
				Ps = np.append(Ps, daug_p)
	# Remove decayed particles
	pdgs = np.delete(pdgs, decayed)
	Es = np.delete(Es, decayed)
	Ps = np.delete(Ps, decayed)
	Pvs = np.delete(Pvs, decayed, 0)

	if verbose:
		print('Event particle IDs', pdgs)
		print('Particle energies', Es)
		print('Particle momenta', Ps)
	# 	print('Particle 3-momenta', Pvs)

	
# True Ring constructor, baseline for reconstructing the event
	TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir = TrueRingConstructor(pdgs, Es, Ps, Pvs)

	if TrueNRing == 0:
		if verbose:
			print('Zero Rings')
		continue


# Reconstructing rings
	# start_time = time.time()

	RRing = RecoRing(TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir, rd, event.Mode[i])

	if RRing.NRing == 1: RRing.DecayE(event.Ipnu[i], event.CC[i], event.Mode[i])

	RRing.SKType(event.Ipnu[i], event.CC[i])

	# print('Recos Ring:', time.time() - start_time, 'seconds')


	if RRing.NRing < 1 or RRing.Type < 0:
		del RRing
		continue

	if verbose:
	# 	if RRing.Nring == 1:
	# 		print('Single Ring event')
	# 		print('Visible energy is ', RRing.Evis, ' GeV')
	# 		print('Event IP is ', RRing.IP)
	# 		print('Number of Michel electrons: ', RRing.MuEdk)
		if RRing.NRing > 1:
			print('Multi Ring event')
			print('Event with ', RRing.NRing, ' reco rings')
			print('Visible energy is ', RRing.Evis, ' GeV')
			print('Event IP is ', RRing.IP)
			print('Event MER IP is ', RRing.MERIP)
			print('Event type is ', RRing.Type)
		# else:
	# 		print('Event with ', RRing.RecoNring, ' reco rings')
	# 		print('Event with ', RRing.TrueNring, ' true rings')

	# Event verbosity
	if verbose:
		print('Neutrino energy: ',event.Enu[i],' GeV')
		print('Neutrino flavour: ',event.Ipnu[i])
		print('Lepton ID: ',event.PDGlep)
		# print('Lepton momentum: ',event.LeptonMomentum, 'GeV')
		if abs(event.Mode[i])<30:
			print('Interacting CC with code: ', event.Mode[i])
		else:
			print('Interacting NC with code: ', event.Mode[i])

	# Fine tuning
	RRing.mendSGE()
	RRing.mendSGM(event.Ipnu[i])
	RRing.mendMG(event.Ipnu[i])
	RRing.mendMR(event.Ipnu[i])

	if RRing.Type<0:
		continue


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
	oscw       = np.append(oscw, event.oscw[i])
	plep       = np.append(plep, event.Plep[i])
	dirlep_x   = np.append(dirlep_x, event.Pxlep[i] / event.Plep[i])
	dirlep_y   = np.append(dirlep_y, event.Pylep[i] / event.Plep[i])
	dirlep_z   = np.append(dirlep_z, event.Pzlep[i] / event.Plep[i])
	# reco variables
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
	mode       = np.append(mode, event.Mode[i])


	del RRing

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
	hf.create_dataset('oscw', data=oscw, compression='gzip')
	hf.create_dataset('itype', data=itype, compression='gzip')
	hf.create_dataset('mode', data=mode, compression='gzip')

