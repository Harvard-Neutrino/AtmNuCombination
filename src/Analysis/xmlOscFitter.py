import numpy as np
import math
from math import asin, sqrt
from SimReader import Reader, OscPar
import multiprocessing
import sys
from XMLreader import parseXML
from itertools import product

# Setup experiment details
analysis_xml_file=str(sys.argv[1])
outfile=str(sys.argv[2])

an = parseXML(analysis_xml_file)
an.readSources()
an.readExperiments()
an.readSystematics()
an.readPhysics()
an.readOscPar()

mc = Reader(an.experiments[0],an.mcFiles[0])
mc.Binning()
mc.SetOutFile(outfile)

# Get unoscillated atm. fluxes
mc.InitialFlux()

# Setup oscillation parameters grid and best fit value
mc.BFOscillatorxml(an.neutrinos,**an.OscParametersBest)

with open(outfile,'w') as f:
	for par in an.parameters:
		f.write(par+' ')
	f.write('X2\n')

mcmc=0

if mcmc==0:
	if an.physics[0] == 'Three Flavour':
		print(an.physics)
		processes = []
		jj = 0

		# Osc. in parameters space (multiprocessing)
		for i,t12 in enumerate(an.OscParametersGrid['Sin2Theta12']):
			for j,t13 in enumerate(an.OscParametersGrid['Sin2Theta13']):
				for k,t23 in enumerate(an.OscParametersGrid['Sin2Theta23']):
					for l,dm21 in enumerate(an.OscParametersGrid['Dm221']):
						for m,dm31 in enumerate(an.OscParametersGrid['Dm231']):
							for d,cp in enumerate(an.OscParametersGrid['dCP']):
								for o,hi in enumerate(an.OscParametersGrid['Ordering']):
									p = multiprocessing.Process(target=mc.Chi2Calculator,args=[an.neutrinos,t12, t13, t23, dm21, dm31, cp, hi])
									if __name__ == "__main__":
										processes.append(p)
										p.start()
										jj = jj + 1
										print(f'{t12} {t13} {t23} {dm21} {dm31} {cp} process started')
										if jj%20==0:
											for i,p in enumerate(processes):
												p.join()
											processes = []
else:
	param = []
	for par in an.parameters:
		if an.OscParametersGrid[par].size>2:
			param.append(par)

	ndim = len(param)
	nwalkers = 32
	print(f'Parameter space of {ndim} dimensions.')
	print(f'We have 32 walkers by default and you are using {nwalkers}.')
	
	p0 = np.random.randn(nwalkers, ndim)

	prior = np.random.rand(ndim)
	prior = an.OscParametersBest

	print(prior)
	for i,par in enumerate(param):
		r = np.random.rand(1)[0]
		prior[par] = an.OscParametersEdges[par][0] + r * (an.OscParametersEdges[par][1] - an.OscParametersEdges[par][0]) 
		# print(an.OscParametersEdges[par])
		# prior[i] = an.OscParametersEdges[par][0] + prior[i] * (an.OscParametersEdges[par][1] - an.OscParametersEdges[par][0]) 

	print(prior)


