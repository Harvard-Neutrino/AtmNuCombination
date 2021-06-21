# Requirements
import ROOT as root
from ROOT import TCanvas, TFile, TTree
import numpy as np
import argparse
#import multiprocessing

# Local dependencies
from RootFilesManager import *
from TrueRing import *
from RecoRing import *
from RandomGenerator import RecoDists
from GenieHDF5 import *
from applications import Azimuth

parser = argparse.ArgumentParser()
parser.add_argument("in_hdf5filename", type=str, nargs='?', default='/home/pablofer/AtmNuCombination/utils/atm_genie.root.hdf5')
parser.add_argument("out_rootfilename", type=str, nargs='?', default='testfcmc.root')
parser.add_argument("-v", dest='verbo', default=False, action='store_true')

args = parser.parse_args()
verbose = args.verbo
genie_input = args.in_hdf5filename
output = args.out_rootfilename

honda_file = '/home/pablofer/IceCube-SuperK_AtmNu/SuperK/Flux/flux_at_Kamioka/honda_flux_log.root'
#distros = 'RecoDistributions.root'

root.gROOT.SetBatch(True) # No graphics enabled

# Initialize random generator
root.gRandom.SetSeed(1)
rand = root.TRandom3(0)

# Read GENIE GST file
event = GenieSimulation('/home/pablofer/AtmNuCombination/utils/atm_genie.root.hdf5')

# Read Honda flux
hfile,htree = read_tree(honda_file,'honda')

# Read reco distributions
#dfile = read_histo(distros)

rd = RecoDists()

# # Output file to write
# ofile, osc_tuple = out_tree(output,'osc_tuple')
# # Variable declaration for out tree
# ipnu       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("ipnu",ipnu,"ipnu/D")
# pnu        = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("pnu",pnu,"pnu/D")
# dirnu_x    = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirnu_x",dirnu_x,"dirnu_x/D")
# dirnu_y    = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirnu_y",dirnu_y,"dirnu_y/D")
# dirnu_z    = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirnu_z",dirnu_z,"dirnu_z/D")
# cz         = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("cz",cz,"cz/D")
# azi        = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("azi",azi,"azi/D")
# fluxho     = np.zeros(4,dtype=np.double) ; osc_tuple.Branch("fluxho",fluxho,"fluxho[4]/D")
# plep       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("plep",plep,"plep/D")
# dirlep_x   = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirlep_x",dirlep_x,"dirlep_x/D")
# dirlep_y   = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirlep_y",dirlep_y,"dirlep_y/D")
# dirlep_z   = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirlep_z",dirlep_z,"dirlep_z/D")
# ## reco variables
# reco_pmax  = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_pmax",reco_pmax,"reco_pmax/D")
# evis       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("evis",evis,"evis/D")
# reco_dir_x = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_dir_x",reco_dir_x,"reco_dir_x/D")
# reco_dir_y = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_dir_y",reco_dir_y,"reco_dir_y/D")
# reco_dir_z = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_dir_z",reco_dir_z,"reco_dir_z/D")
# ip         = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("ip",ip,"ip/D")
# nring      = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("nring",nring,"nring/D")
# muedk      = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("muedk",muedk,"muedk/D")
# itype      = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("itype",itype,"itype/D")
# mode       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("mode",mode,"mode/D")

# if verbose: osc_tuple.Print()

# Start loop over Genie simulated
for i, nu in enumerate(event.Ipnu):

	if i>1000: continue
	
	if verbose:
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

	
# True Ring constructor, baseline for reconstructing the event
	TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir = TrueRingConstructor(event.PDGlep[i], pdg_fp, event.Elep[i], E_fp, event.Plep[i], P_fp, Plep_v, P_fp_v)
	if TrueNRing == 0:
		if verbose:
			print('Zero Rings')
		continue


# Reconstructing rings
	# print(TrueNRing)
	# print(TrueRingPDG)
	# print(TrueRingIP)
	# print(TrueRingE)

	RRing = RecoRing(TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir, rd, event.Mode[i])

	if RRing.NRing == 1: RRing.DecayE(event.Ipnu[i], event.CC[i], event.Mode[i])

	RRing.SKType(event.Ipnu[i], event.CC[i])

	if RRing.NRing < 1 or RRing.Type < 0:
		del RRing
		continue

	# if verbose:
	# 	if RRing.RecoNring == 1:
	# 		print('Single Ring event')
	# 		print('Visible energy is ', RRing.Evis, ' GeV')
	# 		print('Event IP is ', RRing.IP)
	# 		print('Number of Michel electrons: ', RRing.MuEdk)
	# 	elif RRing.RecoNring > 1:
	# 		print('Multi Ring event')
	# 		print('Event with ', RRing.RecoNring, ' reco rings')
	# 		print('Event with ', RRing.TrueNring, ' true rings')
	# 		print('Visible energy is ', RRing.Evis, ' GeV')
	# 		print('Event IP is ', RRing.IP)
	# 	else:
	# 		print('Event with ', RRing.RecoNring, ' reco rings')
	# 		print('Event with ', RRing.TrueNring, ' true rings')

	# # Event verbosity
	# if verbose:
	# 	print('Neutrino energy: ',nu.Energy,' GeV')
	# 	print('Neutrino flavour: ',nu.Flavour)
	# 	print('Lepton ID: ',event.LeptonPDG)
	# 	print('Lepton momentum: ',event.LeptonMomentum, 'GeV')
	# 	if abs(nu.Mode)<30:
	# 		print('Interacting CC with code: ', nu.Mode)
	# 	else:
	# 		print('Interacting NC with code: ', nu.Mode)

	# mp_flux.join()

	ipnu[0]       = event.Ipnu[i]
	pnu[0]        = event.Enu[i]
	dirnu_x[0]    = Pnu_v[0] / event.Enu[i]
	dirnu_y[0]    = Pnu_v[1] / event.Enu[i]
	dirnu_z[0]    = Pnu_v[2] / event.Enu[i]
	cz[0]         = dirnu_z[0]
	azi[0]        = Azimuth(dirnu_y[0], dirnu_z[0])
	fluxho        = event.Flux
	plep[0]       = event.Plep[i]
	dirlep_x[0]   = event.Pxlep[i] / event.Plep[i]
	dirlep_y[0]   = event.Pylep[i] / event.Plep[i]
	dirlep_z[0]   = event.Pzlep[i] / event.Plep[i]
	## reco variables
	reco_pmax[0]  = RRing.MERMomentum
	evis[0]       = RRing.Evis
	reco_dir_x[0] = RRing.TotDir[0]
	reco_dir_y[0] = RRing.TotDir[1]
	reco_dir_z[0] = RRing.TotDir[2]
	ip[0]         = RRing.MERIP
	nring[0]      = RRing.NRing
	muedk[0]      = RRing.MuEdk
	itype[0]      = RRing.Type
	mode[0]       = event.Mode[i]

	del RRing

	osc_tuple.Fill()

osc_tuple.Write()
ofile.Write()
ofile.Close()
