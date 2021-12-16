import numpy as np
import math
from math import asin, sqrt

	
def Chi2StatsCombined(neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile):
	X2 = 0
	for exp in experiments:
		wOsc = experiments[exp].Oscillator(neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering)
		# print(exp)
		for s in range(experiments[exp].NumberOfSamples):
			# print(s)
			wBF = experiments[exp].weightOscBF[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			w = wOsc[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			# print(w-wBF)
			Exp, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=wBF)
			Obs, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=w)
			for O,E in zip(np.ravel(Obs),np.ravel(Exp)):
				# if O==0 or E==0: 
					# print(f'Warning: null statistics at {exp} sample {s}, {E} {O}')
					# pass
				if O!=0 and E!=0:
					X2 = X2 + 2 * (E - O + O * math.log(O/E))
					# if abs(E-O):
					# 	print(f'Warning: different statistics at {exp} sample {s}, {E} {O}')
					# print(X2)
	# print(X2)
	with open(outfile,'a') as f:
		f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {X2}\n')
		f.flush()


	return(X2)


def Chi2SystsCombined(syst, syst_key, neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile):
	X2 = 0
	for exp in experiments:
		wOsc = experiments[exp].Oscillator(neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering)
		for sys in syst:


		for s in range(experiments[exp].NumberOfSamples):
			# print(s)
			wBF = experiments[exp].weightOscBF[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			w = wOsc[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			# print(w-wBF)
			Exp, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=wBF)
			Obs, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=w)
			for O,E in zip(np.ravel(Obs),np.ravel(Exp)):
				# if O==0 or E==0: 
					# print(f'Warning: null statistics at {exp} sample {s}, {E} {O}')
					# pass
				if O!=0 and E!=0:
					X2 = X2 + 2 * (E - O + O * math.log(O/E))
					# if abs(E-O):
					# 	print(f'Warning: different statistics at {exp} sample {s}, {E} {O}')
					# print(X2)
	# print(X2)
	with open(outfile,'a') as f:
		f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {X2}\n')
		f.flush()


	return(X2)
