from scipy.optimize import minimize
from merger import *

def sensitivity(analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile):
	if analysis.NoSyst:
		print('Stats. only')
		Chi2StatsCombined(analysis.neutrinos, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile)

	else:
		bnds = []
		for source in analysis.Systematics:
			for i,sys in enumerate(analysis.Systematics[source]):
				if abs(analysis.SystNominal[source][i] -1) < 0.1:
					bnds.append((0.2,2))
				elif abs(analysis.SystNominal[source][i]) < 0.1:
					bnds.append((-1,1))
		bounds = tuple(bnds)
		res = minimize(Chi2SystsCombined, analysis.SystNominalList, args=(analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments), method='L-BFGS-B', bounds=bounds, tol=0.001)

		# Writing output
		sys_data = ' '.join(map(str,res.x))
		with open(outfile,'a') as f:
			f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {sys_data} {res.fun}\n')
			f.flush()