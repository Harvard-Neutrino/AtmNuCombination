from scipy.optimize import minimize
from merger import Chi2StatsCombined, Chi2SystsCombined, AnalyticPriorsBounds, Chi2StatsBinned, Chi2SystsBinned # , AnalyticPriors
import numpy as np

# def sensitivity(t12, t13, t23, dm21, dm31, dcp, Ordering, analysis, experiments, outfile):
def sensitivity(x, Ordering, analysis, experiments, outfile):
	t12, t13, t23, dm21, dm31, dcp = x
	print(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} process started')
	if t23>1 or t23<0:
		return -9999.

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
		return - 0.5 * statX2

	else:
		# Analytic estimate for priors and bounds
		analysis.SystPrior, bounds = AnalyticPriorsBounds(analysis, Obs, experiments)
		# AnalyticPriors(analysis, Obs, experiments)

		# Combined chi^2 minimization
		tol = max(1e-4,np.sqrt(statX2)*1e-4)
		res = minimize(Chi2SystsCombined, analysis.SystPrior, args=(analysis, Obs, experiments), method='L-BFGS-B', jac=True, bounds=bounds, options={'disp' : False, 'ftol' : tol, 'gtol': 1e-03})

		# Writing output
		sys_data = ' '.join(map(str,res.x))
		with open(outfile,'a') as f:
			f.write(f'{t12:.3f} {t13:.3f} {t23:.5f} {dm21:.7f} {dm31:.7f} {dcp:.3f} {Ordering} {sys_data} {res.fun:.2f}\n')
			f.flush()

		# print(res.fun)
		
		return - 0.5 * res.fun



def binned_sensitivity(x, Ordering, analysis, experiments, outfile):
	t12, t13, t23, dm21, dm31, dcp = x
	print(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} process started')
	if t23>1 or t23<0:
		return -9999.

	# Compute expectation at BF point
	Obs = {}
	for exp in experiments:
		Obs[exp] = experiments[exp].BinOscillator(analysis.neutrinos, t12, t13, t23, dm21, dm31, dcp, Ordering)

	statX2 = Chi2StatsCombined(analysis, Obs, experiments)

	binned_statX2 = Chi2StatsBinned(analysis, Obs, experiments)

	if analysis.NoSyst:
		# Writing output
		with open(outfile,'w') as f:
			counter = -1
			for exp in experiments:
				f.write('Sample E_low E_up Cz_low Cz_up X2\n')
				for sample in range(experiments[exp].NumberOfSamples):
					for energy_bins in range(len(experiments[exp].EnergyBins[sample])-1):
						for cz_bins in range(len(experiments[exp].CzBins[sample])-1):
							counter += 1
							f.write(f'{sample} {experiments[exp].EnergyBins[sample][energy_bins]} {experiments[exp].EnergyBins[sample][energy_bins+1]} {experiments[exp].CzBins[sample][cz_bins]} {experiments[exp].CzBins[sample][cz_bins+1]} {binned_statX2[counter]}')
							f.flush()

		# print(statX2)
		return - 0.5 * statX2

	else:
		# Analytic estimate for priors and bounds
		analysis.SystPrior, bounds = AnalyticPriorsBounds(analysis, Obs, experiments)
		# AnalyticPriors(analysis, Obs, experiments)

		# Combined chi^2 minimization
		tol = max(1e-4,np.sqrt(statX2)*1e-4)
		res = minimize(Chi2SystsCombined, analysis.SystPrior, args=(analysis, Obs, experiments), method='L-BFGS-B', jac=True, bounds=bounds, options={'disp' : False, 'ftol' : tol, 'gtol': 1e-03})

		# Writing output
		sys_data = ' '.join(map(str,res.x))
		# with open(outfile,'a') as f:
		# 	f.write(f'{t12:.3f} {t13:.3f} {t23:.5f} {dm21:.7f} {dm31:.7f} {dcp:.3f} {Ordering} {sys_data} {res.fun:.2f}\n')
		# 	f.flush()
		binned_systX2 = Chi2SystsBinned(res.x, analysis, Obs, experiments)
		with open(outfile,'w') as f:
			counter = -1
			for exp in experiments:
				f.write('Sample E_low E_up Cz_low Cz_up X2\n')
				for sample in range(experiments[exp].NumberOfSamples):
					for energy_bins in range(len(experiments[exp].EnergyBins[sample])-1):
						for cz_bins in range(len(experiments[exp].CzBins[sample])-1):
							counter += 1
							f.write(f'{sample} {experiments[exp].EnergyBins[sample][energy_bins]} {experiments[exp].EnergyBins[sample][energy_bins+1]} {experiments[exp].CzBins[sample][cz_bins]} {experiments[exp].CzBins[sample][cz_bins+1]} {binned_systs2[counter]}')
							f.flush()

		# print(res.fun)
		
		return - 0.5 * res.fun

