# Requirements
import ROOT as root
from ROOT import TCanvas, TFile, TTree
import numpy as np
import math
import argparse
#import multiprocessing

# Local dependencies
from RootFilesManager import *
from Event import *
from Neutrino import *
from TrueRing import *
from RecoRing import *
from applications import ClearBadHadrons
from RandomGenerator import RecoDists

parser = argparse.ArgumentParser()
parser.add_argument("in_rootfilename", type=str, nargs='?', default='/home/pablofer/IceCube-SuperK_AtmNu/SuperK/AtmNu_MC/atm_genie.root')
parser.add_argument("out_rootfilename", type=str, nargs='?', default='fcmc.root')
parser.add_argument("-v", dest='verbo', default=False, action='store_true')

args = parser.parse_args()
verbose = args.verbo
genie_input = args.in_rootfilename
output = args.out_rootfilename

honda_file = '/home/pablofer/IceCube-SuperK_AtmNu/SuperK/Flux/flux_at_Kamioka/honda_flux_log.root'
#distros = 'RecoDistributions.root'

root.gROOT.SetBatch(True) # No graphics enabled

# Initialize random generator
root.gRandom.SetSeed(1)
rand = root.TRandom3(0)

# Read GENIE GST file
gfile, gtree = read_tree(genie_input,'gst')

# Read Honda flux
hfile,htree = read_tree(honda_file,'honda')

# Read reco distributions
#dfile = read_histo(distros)

rd = RecoDists()

# Output file to write
ofile, osc_tuple = out_tree(output,'osc_tuple')
# Variable declaration for out tree
ipnu       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("ipnu",ipnu,"ipnu/D")
pnu        = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("pnu",pnu,"pnu/D")
dirnu_x    = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirnu_x",dirnu_x,"dirnu_x/D")
dirnu_y    = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirnu_y",dirnu_y,"dirnu_y/D")
dirnu_z    = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirnu_z",dirnu_z,"dirnu_z/D")
cz         = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("cz",cz,"cz/D")
azi        = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("azi",azi,"azi/D")
fluxho     = np.zeros(4,dtype=np.double) ; osc_tuple.Branch("fluxho",fluxho,"fluxho[4]/D")
plep       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("plep",plep,"plep/D")
dirlep_x   = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirlep_x",dirlep_x,"dirlep_x/D")
dirlep_y   = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirlep_y",dirlep_y,"dirlep_y/D")
dirlep_z   = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("dirlep_z",dirlep_z,"dirlep_z/D")
## reco variables
reco_pmax  = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_pmax",reco_pmax,"reco_pmax/D")
evis       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("evis",evis,"evis/D")
reco_dir_x = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_dir_x",reco_dir_x,"reco_dir_x/D")
reco_dir_y = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_dir_y",reco_dir_y,"reco_dir_y/D")
reco_dir_z = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("reco_dir_z",reco_dir_z,"reco_dir_z/D")
ip         = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("ip",ip,"ip/D")
nring      = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("nring",nring,"nring/D")
muedk      = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("muedk",muedk,"muedk/D")
itype      = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("itype",itype,"itype/D")
mode       = np.zeros(1,dtype=np.double) ; osc_tuple.Branch("mode",mode,"mode/D")

if verbose: osc_tuple.Print()

if verbose: print('Analyzing ',gtree.GetEntries(),' SK atmospheric neutrinos,')

# Start loop over Genie simulated
for i, entry in enumerate(gtree):

	# if i>1000: continue
	
	if verbose:
	    print('----------------------------------------------------')
	    print('----------- Event number ',i,'--------------------')
	    print('----------------------------------------------------')

	# Get variables from Genie tree
	# Neutrino
	Pnu_v=np.array([entry.pxv,entry.pyv,entry.pzv]) # Neutrino momentum vector
	# Primary lepton
	Plep_v=np.array([entry.pxl,entry.pyl,entry.pzl]) # Prompt lepton momentum vector

	# Hadronic system arrangement
	n_fp = entry.nf
	pdg_fp=(np.array(entry.pdgf)) # PDG IDs
	E_fp=(np.array(entry.Ef)) # Energies
	dP_fp_v=np.array([np.array(entry.pxf),np.array(entry.pyf),np.array(entry.pzf)]) # Momenta
	P_fp_v = np.transpose(dP_fp_v)
	P_fp=np.transpose(np.linalg.norm(dP_fp_v, axis=0)) # Momentum
	Dirfp = np.zeros((n_fp,3))
	for j in range(n_fp):
		for k in range(3):
			Dirfp[j][k] = P_fp_v[j][k] / P_fp[j]
	norm = np.linalg.norm(Dirfp, axis=1)

	# Clean up hadrons with null momentum
	n_fp, norm, E_fp, pdg_fp, Dirfp, P_fp_v, P_fp = ClearBadHadrons(n_fp, norm, E_fp, pdg_fp, Dirfp, P_fp_v, P_fp)

# Neutrino class
	nu = Neutrino(htree, entry.neu, entry.Ev, Pnu_v, entry.neut_code, entry.cc, entry.coh, entry.dis, entry.nc)
	# Compute fluxes for this neutrino and its flavour companions
	HKKMflux = nu.Flux()
	# mp_flux = multiprocessing.Process(target=Flux, args=[nu])
	# mp_flux.start()
	# print(HKKMflux)
	
# True event class
	event = Event(nu.LeptonPDG(),entry.El,Plep_v,entry.pl,entry.nf,pdg_fp,E_fp,Dirfp,P_fp_v,P_fp)
	if event.NumberOfParticles == 0:
		if verbose:
			print('Zero particles: ??')
		continue

# True rings
	true_ring = TrueRing(event, nu)
	if true_ring.Nring < 1:
		if verbose:
			print('Zero true rings')
		continue
		

# Reconstructing rings
	#reco_ring = RecoRing(true_ring, dfile)
	reco_ring = RecoRing(true_ring, rd)
	if reco_ring.RecoNring == 1: reco_ring.DecayE(nu.Flavour, nu.Mode)
	reco_ring.SKType(nu.Flavour, nu.Mode)

	if reco_ring.RecoNring < 1:
		continue

	if reco_ring.Type < 0:
		continue

	if verbose:
		if reco_ring.RecoNring == 1:
			print('Single Ring event')
			print('Visible energy is ', reco_ring.Evis, ' GeV')
			print('Event IP is ', reco_ring.IP)
			print('Number of Michel electrons: ', reco_ring.MuEdk)
		elif reco_ring.RecoNring > 1:
			print('Multi Ring event')
			print('Event with ', reco_ring.RecoNring, ' reco rings')
			print('Event with ', reco_ring.TrueNring, ' true rings')
			print('Visible energy is ', reco_ring.Evis, ' GeV')
			print('Event IP is ', reco_ring.IP)
		else:
			print('Event with ', reco_ring.RecoNring, ' reco rings')
			print('Event with ', reco_ring.TrueNring, ' true rings')

	# Event verbosity
	if verbose:
		print('Neutrino energy: ',nu.Energy,' GeV')
		print('Neutrino flavour: ',nu.Flavour)
		print('Lepton ID: ',event.LeptonPDG)
		print('Lepton momentum: ',event.LeptonMomentum, 'GeV')
		if abs(nu.Mode)<30:
			print('Interacting CC with code: ', nu.Mode)
		else:
			print('Interacting NC with code: ', nu.Mode)

	# mp_flux.join()
	ipnu[0]       = nu.Flavour
	pnu[0]        = nu.Energy
	dirnu_x[0]    = nu.Direction[0]
	dirnu_y[0]    = nu.Direction[1]
	dirnu_z[0]    = nu.Direction[2]
	cz[0]         = nu.CosZenith
	azi[0]        = nu.Azi
	fluxho        = HKKMflux
	plep[0]       = event.LeptonMomentum
	dirlep_x[0]   = event.LeptonDirection[0]
	dirlep_y[0]   = event.LeptonDirection[1]
	dirlep_z[0]   = event.LeptonDirection[2]
	## reco variables
	reco_pmax[0]  = reco_ring.MERMomentum
	evis[0]       = reco_ring.Evis
	reco_dir_x[0] = reco_ring.Direction[0]
	reco_dir_y[0] = reco_ring.Direction[1]
	reco_dir_z[0] = reco_ring.Direction[2]
	ip[0]         = reco_ring.IP
	nring[0]      = reco_ring.RecoNring
	muedk[0]      = reco_ring.MuEdk
	itype[0]      = reco_ring.Type
	mode[0]       = nu.Mode

	osc_tuple.Fill()

osc_tuple.Write()
ofile.Write()
ofile.Close()
