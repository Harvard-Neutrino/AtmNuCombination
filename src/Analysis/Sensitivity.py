from scipy.optimize import minimize
from merger import *
import numpy as np

def sensitivity(analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile):
	
	if analysis.NoSyst:
		Chi2StatsCombined(analysis.neutrinos, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile)

	else:
		bnds = []
		prior = []
		for source in analysis.Systematics:
			for i,sys in enumerate(analysis.Systematics[source]):
				# prior.append(analysis.SystNominal[source][i]+0.02*np.random.rand())
				prior.append(analysis.SystNominal[source][i])
				if abs(analysis.SystNominal[source][i] -1) < 0.1:
					bnds.append((0.5,2))
				elif abs(analysis.SystNominal[source][i]) < 0.1:
					bnds.append((-0.9,0.9))
		bounds = tuple(bnds)
		res = minimize(Chi2SystsCombined, prior, args=(analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments), method='L-BFGS-B', bounds=bounds, tol=0.001)

		# Writing output
		sys_data = ' '.join(map(str,res.x))
		with open(outfile,'a') as f:
			f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {sys_data} {res.fun}\n')
			f.flush()