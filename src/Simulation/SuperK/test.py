# Requirements
# import ROOT as root
# from ROOT import TCanvas, TFile, TTree
import numpy as np
import math
import argparse
#import multiprocessing

# Local dependencies
# from RootFilesManager import *
# from Event import *
# from Neutrino import *
# from TrueRing import *
# from RecoRing import *
# from applications import ClearBadHadrons
# from RandomGenerator import RecoDists
from GenieHDF5 import *

event = GenieSimulation('/home/pablofer/AtmNuCombination/utils/atm_genie.root.hdf5')
for i, nu in enumerate(event.Ipnu):
	# indeces = event.Phad[i] < 0.001
	# P_fp_v = np.transpose(np.array([event.Pxhad[i], event.Pyhad[i], event.Pzhad[i]])) # Momenta
	slowHad = event.RestHadrons(i)
	if not np.any(slowHad):
		px = event.Pxhad[i]
		py = event.Pyhad[i]
		pz = event.Pzhad[i]
		p = event.Phad[i]
		e = event.Ehad[i]
		pdg = event.PDGhad[i] # PDG IDs
		P_fp_v = np.array([px[slowHad], py[slowHad], pz[slowHad]]) # Momenta
		P_fp = p[slowHad]
		E_fp = e[slowHad]
		pdg_fp = pdg[slowHad]
		n_fp = np.size(P_fp)
	else:
		P_fp_v = np.array([event.Pxhad[i], event.Pyhad[i], event.Pzhad[i]]) # Momenta
		P_fp = event.Phad[i]
		n_fp = event.Nhad
		pdg_fp = event.PDGhad[i] # PDG IDs
		E_fp = event.Ehad[i] # Energies
		

	Dirfp = P_fp_v / P_fp # Direction
	P_fp_v = np.transpose(P_fp_v)
	Dirfp = np.transpose(Dirfp)
	print(Dirfp)
	print('---------------------------------')

