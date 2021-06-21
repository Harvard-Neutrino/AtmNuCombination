# Requirements
import ROOT as root
from ROOT import TCanvas, TFile, TTree
import numpy as np
import math
import argparse
#import multiprocessing

# Local dependencies
from RootFilesManager import *

parser = argparse.ArgumentParser()
parser.add_argument("in_rootfilename", type=str, nargs='?', default='/home/pablofer/IceCube-SuperK_AtmNu/SuperK/AtmNu_MC/atm_genie.root')

args = parser.parse_args()
genie_input = args.in_rootfilename

root.gROOT.SetBatch(True) # No graphics enabled

# Read GENIE GST file
gfile, gtree = read_tree(genie_input,'gst')
for i, entry in enumerate(gtree):
	j=entry.iev