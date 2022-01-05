import numpy as np
import math
from math import asin, sqrt
from SimReader import Reader
import multiprocessing
import sys
from xmlReader import parseXML
from itertools import product
from merger import Chi2StatsCombined
from Sensitivity import sensitivity

# Setup analysis details
analysis_xml_file=str(sys.argv[1])
outfile=str(sys.argv[2])

an = parseXML(analysis_xml_file)
an.readSources()
an.readExperiments()
an.readPhysics()
an.readOscPar()
an.CheckSystematics()

# Analysis flags
mcmc = 0
multiproc = 1

# Setup all experiments
mcList = {}

for s in an.sources:
	for i,(exp,fil) in enumerate(zip(an.experiments,an.mcFiles)):
		mcList[exp] = Reader(s,exp,fil)
		mcList[exp].Binning()
		# Get unoscillated atm. fluxes
		mcList[exp].InitialFlux()
		# Setup oscillation parameters grid and best fit value
		mcList[exp].BFOscillator(an.neutrinos,**an.OscParametersBest)


print('=============================================================')
print('====================== There we go!! ========================')


with open(outfile,'w') as f:
	for par in an.parameters:
		f.write(par+' ')
		if an.NoSyst == 0:
			for sys in an.Systematics.values():
				for s in sys:
					f.write(s+' ')
		f.write('X2 ')
		f.write('\n')


# Main analysis

if multiproc:
	cores = multiprocessing.cpu_count()

	if mcmc==0:
		if an.physics[0] == 'Three Flavour':
			processes = []
			jj = 0

			# Osc. in parameters space (multiprocessing)
			param = []
			for oscpar in [*an.OscParametersGrid]:
				param.append(an.OscParametersGrid[oscpar])

			for element in product(*param):
				p = multiprocessing.Process(target=sensitivity,args=[an, *element, mcList, outfile])
				if __name__ == "__main__":
					processes.append(p)
					p.start()
					jj = jj + 1
					print(f'{element} process started')
					if jj%cores==0:
						for i,p in enumerate(processes):
							p.join()
							processes = []
							print('------------------------------')

	'''
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


	'''
else:
	if mcmc==0:
		if an.physics[0] == 'Three Flavour':
		
			# Osc. in parameters space (multiprocessing)
			param = []
			for oscpar in [*an.OscParametersGrid]:
				param.append(an.OscParametersGrid[oscpar])

			for element in product(*param):
				sensitivity(an, *element, mcList, outfile)

	'''
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


	'''
