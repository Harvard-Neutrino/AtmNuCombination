from scipy.optimize import minimize
from merger import Chi2StatsCombined, Chi2SystsCombined, AnalyticPriorsBounds
import numpy as np

def sensitivity(analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile):
	
	# Compute expectation at BF point
	Obs = {}
	for exp in experiments:
		Obs[exp] = experiments[exp].BinOscillator(analysis.neutrinos, t12, t13, t23, dm21, dm31, dcp, Ordering)

	statX2 = Chi2StatsCombined(analysis, Obs, experiments)

	if analysis.NoSyst:
		# Writing output
		with open(outfile,'a') as f:
			f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {statX2}\n')
			f.flush()

		# print(statX2)
		return statX2

	else:
		# Analytic estimate for priors and bounds
		analysis.SystPrior, bounds = AnalyticPriorsBounds(analysis, Obs, experiments)

		# Combined chi^2 minimization
		tol = max(1e-4,np.sqrt(statX2)*1e-4)
		res = minimize(Chi2SystsCombined, analysis.SystPrior, args=(analysis, Obs, experiments), method='L-BFGS-B', jac=True, bounds=bounds, options={'disp' : False, 'ftol' : tol, 'gtol': 1e-03})

		# Writing output
		sys_data = ' '.join(map(str,res.x))
		with open(outfile,'a') as f:
			f.write(f'{t12:.3f} {t13:.3f} {t23:.5f} {dm21:.7f} {dm31:.7f} {dcp:.3f} {Ordering} {sys_data} {res.fun:.2f}\n')
			f.flush()

		# print(res.fun)
		
		return res.fun




