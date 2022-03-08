from scipy.optimize import minimize
from merger import *
import numpy as np

def sensitivity(analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile):
	
	if analysis.NoSyst:
		Chi2StatsCombined(analysis.neutrinos, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile)

	else:
		bnds = []
		prior = []
		for prior,sigma in zip(analysis.SystPrior,analysis.SystSigmaList):
			if sigma<0.05:
				bnd_min = prior - 2*sigma
				bnd_max = prior + 2*sigma
			else:
				bnd_min = prior - sigma
				bnd_max = prior + sigma
			if prior>0:
				bnd_min = max(0.5,bnd_min)
				bnd_max = min(1.5,bnd_max)
			elif prior==0:
				bnd_min = max(-0.5,bnd_min)
				bnd_max = min(0.5,bnd_max)

			bnds.append((bnd_min,bnd_max))
		bounds = tuple(bnds)

		wOsc = {}
		for exp in experiments:
			wOsc[exp] = experiments[exp].Oscillator(analysis.neutrinos, t12, t13, t23, dm21, dm31, dcp, Ordering)

		# res = minimize(Chi2SystsCombined_and_Gradient, prior, args=(analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments), method='L-BFGS-B', jac=True, bounds=bounds, options={'disp' : True, 'ftol' : 1e-4})
		res = minimize(Chi2SystsCombined_and_Gradient_Binned_Next, analysis.SystPrior, args=(analysis, wOsc, experiments), method='L-BFGS-B', jac=True, bounds=bounds, options={'disp' : True, 'ftol':1e-4, 'gtol': 1e-01})

		analysis.SystPrior = res.x

		# Writing output
		sys_data = ' '.join(map(str,res.x))
		with open(outfile,'a') as f:
			f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {sys_data} {res.fun}\n')
			f.flush()