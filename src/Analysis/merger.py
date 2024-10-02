import numpy as np
import math
from scipy.special import gamma
from Systematics import *
	
def Chi2StatsCombined(analysis, Obs, experiments):
	X2 = 0
	for exp in experiments.values():
		# Binned tatistics
		E = exp.weightOscBF_binned
		O = Obs[exp.Experiment]
		X2 += 2 * np.sum(E-O+O*np.log(O/E))

	return X2


def Chi2StatsBinned(analysis, Obs, experiments):
	for exp in experiments.values():
		# Binned tatistics
		E = exp.weightOscBF_binned
		O = Obs[exp.Experiment]
		X2 = 2 * (E-O+O*np.log(O/E))

	return X2


# def AnalyticPriors(analysis, Obs, experiments):
# 	priors = [0] * len(analysis.SystPrior)
# 	A = [0] * len(analysis.SystPrior)
# 	B = [0] * len(analysis.SystPrior)
# 	A = np.zeros(len(analysis.SystPrior))
# 	B = np.zeros((len(analysis.SystPrior),len(analysis.SystPrior)))

# 	# Experiments
# 	for exp in experiments.values():

# 		# Binned statistics
# 		E = exp.weightOscBF_binned
# 		O = Obs[exp.Experiment]
# 		OmE = O-E
# 		EmO = - OmE

# 		#Systematics
# 		usedSysts = []
# 		dfdx = [0] * len(analysis.SystPrior)
# 		thisSyst = analysis.Systematics[exp.Experiment] + analysis.Systematics[exp.Source] + analysis.Systematics[exp.Detector]
# 		for sys in thisSyst:
# 			index = np.where(analysis.SystematicsList==sys)[0]
# 			j = index[0]
# 			dfdx[j] = globals()['Diff_'+sys](analysis.SystNominalList[j],exp)
# 			usedSysts.append(j)

# 		# Compute Jacobian of Chi^2
# 		for i in usedSysts:
# 			A[i] += np.sum(EmO * dfdx[i])
# 			B[i] += np.sum(( dfdx[i] * dfdx)) 
# 			B[i,i] += 1/analysis.SystSigmaList[i]**2

# 	print(np.around(A))
# 	print(np.around(B))
# 	print(analysis.SystNominalList - np.matmul(np.linalg.inv(B),A.T))


def AnalyticPriorsBounds(analysis, Obs, experiments):
	priors = [0] * len(analysis.SystPrior)
	A = [0] * len(analysis.SystPrior)
	B = [0] * len(analysis.SystPrior)

	# Experiments
	for exp in experiments.values():

		# Binned statistics
		E = exp.weightOscBF_binned
		O = Obs[exp.Experiment]
		OmE = O-E

		#Systematics
		usedSysts = []
		dfdx = [0] * len(analysis.SystPrior)
		thisSyst = analysis.Systematics[exp.Experiment] + analysis.Systematics[exp.Source] + analysis.Systematics[exp.Detector]
		for sys in thisSyst:
			index = np.where(analysis.SystematicsList==sys)[0]
			j = index[0]
			dfdx[j] = globals()['Diff_'+sys](analysis.SystNominalList[j],exp)
			usedSysts.append(j)

		# Compute Jacobian of Chi^2
		for i in usedSysts:
			A[i] += np.sum(OmE * dfdx[i])
			B[i] += np.sum(O * dfdx[i] * dfdx[i])

	# Systematic's penalty terms
	bnds = []
	for i,(mu,sig) in enumerate(zip(analysis.SystNominalList,analysis.SystSigmaList)):
		pr = mu + A[i] / (B[i] + 1/sig**2)
		delta = min(np.abs(pr-mu), sig)
		priors[i] = 0.5*(mu+pr)
		if delta>0:
			bnds.append((priors[i]-delta, priors[i]+delta))
		else:
			bnds.append((priors[i]-sig, priors[i]+sig))
	return priors, tuple(bnds)

def Chi2SystsCombined(syst, analysis, Obs, experiments):
	JX2 = [0] * len(syst)
	X2 = 0

	# Experiments
	for exp in experiments.values():

		# Binned tatistics
		E = exp.weightOscBF_binned
		O = Obs[exp.Experiment]

		#Systematics
		usedSysts = []
		dEdx = [0] * len(syst)
		wSys = 0
		dummywSys = 0
		thisSyst = analysis.Systematics[exp.Experiment] + analysis.Systematics[exp.Source] + analysis.Systematics[exp.Detector]
		for sys in thisSyst:
			index = np.where(analysis.SystematicsList==sys)[0]
			j = index[0]
			xFij = globals()[sys](syst[j],exp)
			wSys += xFij
			dEdx[j] = E * globals()['Diff_'+sys](syst[j],exp)
			usedSysts.append(j)

		Es = E * (1 + wSys)

		# Compute Chi^2
		if np.any(Es<=0):
			X2 = 1e6
		else:
			X2 += 2 * np.sum(Es-O+O*np.log(O/Es))

		# Compute Jacobian of Chi^2
		for i in usedSysts:
			JX2[i] += 2 * np.sum((1-O/Es)*dEdx[i])

	# Systematic's penalty terms
	for i,(x,mu,sig) in enumerate(zip(syst,analysis.SystNominalList,analysis.SystSigmaList)):
		X2 += ((x-mu) / sig)**2
		JX2[i] += 2 * (x-mu) / sig**2

	return (X2,JX2)

def Chi2SystsBinned(syst, analysis, Obs, experiments):
	# Experiments
	for exp in experiments.values():

		# Binned tatistics
		E = exp.weightOscBF_binned
		O = Obs[exp.Experiment]

		#Systematics
		usedSysts = []
		wSys = 0
		dummywSys = 0
		thisSyst = analysis.Systematics[exp.Experiment] + analysis.Systematics[exp.Source] + analysis.Systematics[exp.Detector]
		for sys in thisSyst:
			index = np.where(analysis.SystematicsList==sys)[0]
			j = index[0]
			xFij = globals()[sys](syst[j],exp)
			wSys += xFij
			usedSysts.append(j)

		Es = E * (1 + wSys)

		# Compute Chi^2 for each bin
		X2 = 2 * (Es-O+O*np.log(O/Es))

	# Systematic's penalty terms are ignored for this calculation as they are spreaded over all bins
	# for i,(x,mu,sig) in enumerate(zip(syst,analysis.SystNominalList,analysis.SystSigmaList)):
	# 	X2 += ((x-mu) / sig)**2

	return X2
