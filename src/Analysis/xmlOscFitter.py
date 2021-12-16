import numpy as np
import math
from math import asin, sqrt
from SimReader import Reader
import multiprocessing
import sys
from XMLreader import parseXML
from itertools import product
from merger import Chi2StatsCombined

# Setup analysis details
analysis_xml_file=str(sys.argv[1])
if len(sys.argv)==4:
	outfile=[str(sys.argv[2]),str(sys.argv[3])]
elif len(sys.argv)==3:
	outfile=[str(sys.argv[2])]
elif len(sys.argv)==5:
	outfile=[str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4])]
	comboutfile = str(sys.argv[4])
else:
	outfile=str(sys.argv[2])


an = parseXML(analysis_xml_file)
an.readSources()
an.readExperiments()
an.readSystematics()
an.readPhysics()
an.readOscPar()

# Setup all experiments
mcList = {}
# with open(outfile,'w') as f:
# 	for par in an.parameters:
# 		f.write(par+' '+str(an.OscParametersBest[par])+' ')
# 	f.write('X2\n')

# print(len(outfile))

for i,(exp,fil) in enumerate(zip(an.experiments,an.mcFiles)):
	print(exp,fil)
	mcList[exp] = Reader(exp,fil)
	mcList[exp].Binning()
	# mcList[exp].SetOutFile(outfile[0])
	# mcList[exp].SetOutFile(outfile[i])

	# Get unoscillated atm. fluxes
	mcList[exp].InitialFlux()

	# Setup oscillation parameters grid and best fit value
	mcList[exp].BFOscillator(an.neutrinos,**an.OscParametersBest)


# print(an.OscParametersBest)
# chi2Comb = Chi2StatsCombined(3, an.OscParametersBest['Sin2Theta12'], an.OscParametersBest['Sin2Theta13'], an.OscParametersBest['Sin2Theta23'], 
# 		an.OscParametersBest['Dm221'], an.OscParametersBest['Dm231'], an.OscParametersBest['dCP'], an.OscParametersBest['Ordering'], mcList)


if len(outfile)>1:
	for fname in outfile:
		with open(fname,'w') as f:
			for par in an.parameters:
				f.write(par+' ')
			f.write('X2\n')
else:
	with open(outfile[0],'w') as f:
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
									# print(hi)
									# X0 = Chi2StatsCombined(an.neutrinos,t12, t13, t23, dm21, dm31, cp, hi, mcList, outfile[0])
									# print(X0)
									p = multiprocessing.Process(target=Chi2StatsCombined,args=[an.neutrinos,t12, t13, t23, dm21, dm31, cp, hi, mcList, outfile[0]])
									# for mc in mcList:
										# Xf = Xf + mcList[mc].Chi2Calculator(an.neutrinos,t12, t13, t23, dm21, dm31, cp, hi)
										# p = multiprocessing.Process(target=mcList[mc].Chi2Calculator,args=[an.neutrinos,t12, t13, t23, dm21, dm31, cp, hi])
									if __name__ == "__main__":
										processes.append(p)
										p.start()
										jj = jj + 1
										print(f'{t12} {t13} {t23} {dm21} {dm31} {cp} process started')
										if jj%20==0:
											for i,p in enumerate(processes):
												p.join()
											processes = []
											print('Bunch done')
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
