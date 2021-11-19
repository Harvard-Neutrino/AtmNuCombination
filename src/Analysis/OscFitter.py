import numpy as np
import math
import argparse
from SimReader import Reader, OscPar
import multiprocessing
import sys
# from itertools import product

# Setup experiment details
''' ../Simulation/SuperK/data/output/combined.hdf5'''
''' ../Simulation/IceCubeUpgrade/data/neutrino_mc.csv'''
experiment=str(sys.argv[1])
infile=str(sys.argv[2])
outfile=str(sys.argv[3])
mc = Reader(experiment,infile)
mc.Binning()

# Get unoscillated atm. fluxes
mc.InitialFlux()

# Setup oscillation parameters grid and best fit value
# theta_23
NT23 = 24
SqT23_min = 0.305
SqT23_max = 0.705
# dm_31
NDM31 = 24
NDM31 = 0
DM31_max = 3.0e-3
DM31_min = 2.0e-3
# dcp
NDCP = 19
NDCP = 0
DCP_min = 0
DCP_max = 2*math.pi
DM31,T23,SqT23,DCP = OscPar(NDM31,DM31_min,DM31_max,NT23,SqT23_min,SqT23_max,NDCP,DCP_min,DCP_max)
# Best fit values
T23BF = T23[16]
# DM31BF = DM31[12]
# DCPBF = DCP[13]
DM31BF = DM31[0]
DCPBF = DCP[0]

# Start oscillation bussiness
# mc.OscillationSetUp()
mc.BFOscillator(T23BF,DM31BF,DCPBF)

manager = multiprocessing.Manager()
return_dict = manager.dict()

processes = []
j=0

# Osc. in parameters space (multiprocessing)
for t,t23 in enumerate(T23):
	for m,dm in enumerate(DM31):
		for d,cp in enumerate(DCP):
			p = multiprocessing.Process(target=mc.Oscillator,args=[j,t23,dm,cp,return_dict])
			if __name__ == "__main__":
				processes.append(p)
				p.start()
				j = j + 1
for p in processes:
	p.join()


with open(outfile,'w') as f:
	f.write('SqT23 DM31 DCP X2\n')
	j = 0
	for t in range(NT23+1):
		for m in range(NDM31+1):
			for d in range(NDCP+1):
				X2 = 0.0
				weightOsc = return_dict[j]
				j = j + 1
				for s in range(mc.NumberOfSamples):
					wBF = mc.weightOscBF[mc.Sample==s] * mc.Weight[mc.Sample==s]
					Exp, dx, dy = np.histogram2d(mc.EReco[mc.Sample==s], mc.CosZReco[mc.Sample==s], bins=(mc.EnergyBins[s], mc.CzBins[s]), weights=wBF*mc.Norm)
					w = weightOsc[mc.Sample==s] * mc.Weight[mc.Sample==s]
					Obs, dx, dy = np.histogram2d(mc.EReco[mc.Sample==s], mc.CosZReco[mc.Sample==s], bins=(mc.EnergyBins[s], mc.CzBins[s]), weights=w*mc.Norm)
					for O,E in zip(np.ravel(Obs),np.ravel(Exp)):
						if O==0 or E==0: 
							pass
						else:
							X2 = X2 + 2 * (E - O + O * math.log(O/E))
				f.write(f'{SqT23[t]} {DM31[m]} {DCP[d]} {X2}\n')