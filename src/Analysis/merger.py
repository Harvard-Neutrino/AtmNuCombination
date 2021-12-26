import numpy as np
import math
from Systematics import *

	
def Chi2StatsCombined(neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments, outfile):
	X2 = 0
	for exp in experiments:
		wOsc = experiments[exp].Oscillator(neutrino_flavors, t12, t13, t23, dm21, dm31, dcp, Ordering)
		for s in range(experiments[exp].NumberOfSamples):
			wBF = experiments[exp].weightOscBF[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			w = wOsc[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			Exp, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=wBF)
			Obs, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=w)
			for O,E in zip(np.ravel(Obs),np.ravel(Exp)):
				if O!=0 and E!=0:
					X2 = X2 + 2 * (E - O + O * math.log(O/E))
					# print(f'sample, e, o: {s} {O} {E} --> {X2}')

	# Writing output
	with open(outfile,'a') as f:
		f.write(f'{t12} {t13} {t23} {dm21} {dm31} {dcp} {Ordering} {X2}\n')
		f.flush()




def Chi2SystsCombined(syst, analysis, t12, t13, t23, dm21, dm31, dcp, Ordering, experiments):
	X2 = 0
	for exp in experiments:
		wOsc = experiments[exp].Oscillator(analysis.neutrinos, t12, t13, t23, dm21, dm31, dcp, Ordering)

		# Systematics
		wSys = np.ones(experiments[exp].NumberOfEvents)
		for i,source in enumerate(analysis.Systematics):
			if source == exp or source == experiments[exp].Source:
				for sys in analysis.Systematics[source]:
					index = np.where(analysis.SystematicsList==sys)[0]
					j = index[0]
					wSys = wSys * globals()[sys](syst[j],experiments[exp])
					X2 = X2 + ((syst[j]-analysis.SystNominalList[j])/analysis.SystSigmaList[j])**2

		# Statistics
		for s in range(experiments[exp].NumberOfSamples):
			wBF = experiments[exp].weightOscBF[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			w = wSys[experiments[exp].Sample==s] * wOsc[experiments[exp].Sample==s] * experiments[exp].Weight[experiments[exp].Sample==s] * experiments[exp].Norm
			Exp, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=wBF)
			Obs, dx, dy = np.histogram2d(experiments[exp].EReco[experiments[exp].Sample==s], experiments[exp].CosZReco[experiments[exp].Sample==s], bins=(experiments[exp].EnergyBins[s], experiments[exp].CzBins[s]), weights=w)
			for O,E in zip(np.ravel(Obs),np.ravel(Exp)):
				if O!=0 and E!=0:
					X2 = X2 + 2 * (E - O + O * math.log(O/E))

	return(X2)
