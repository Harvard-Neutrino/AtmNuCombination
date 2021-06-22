import h5py
import numpy as np
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

with h5py.File('../Simulation/SuperK/testfcmc.hdf5', 'r') as hf:
	evis = np.array(hf['evis'][()])
	cz = np.array(hf['recodirZ'][()])
	mode = np.array(hf['mode'][()])
	ipnu = np.array(hf['ipnu'][()])
	itype = np.array(hf['itype'][()])



# Livetime normalizations
years=14.58
livetime = years*365.25
total = np.size(itype)
norm = livetime / (total / 8.3)

# Energy bining
sge_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
sgm_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
sgsrpi0ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
sgmrpi0ebins = np.array([0.1, 0.15, 0.25, 0.4, 0.63, 1.33])
mge_ebins = np.array([1.3, 2.5, 5., 10., 20.])
mgm_ebins = np.array([1.3, 3.0, 20.])
mre0_ebins = np.array([1.3, 3.0, 7.0, 20.])
mre1_ebins = np.array([1.3, 2.5, 5., 10., 20.])
mrm_ebins = np.array([1.3, 2.5, 5., 10., 20.])

# Zenith bining
z10bins = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
z1bins = np.array([-1, 1.0])

# Names of SK FC samples
SampleLabel = np.array(['SingleRingSubGeV Elike 0dcy-e','SingleRing SubGeV Elike 1dcy-e','SingleRing Pi0like','SingleRing SubGeV Mulike 0dcy-e',
	'SingleRing SubGeV Mulike 1dcy-e','SingleRing SubGeV Mulike 2dcy-e','MultiRing Pi0like','SingleRing MultiGeV Elike Nulike','SingleRing MultiGeV Elike NuBarlike',
	'SingleRing MultiGeV Mulike','MultiRing MultiGeV Elike Nulike','MultiRing MultiGeV Elike NuBarlike','MultiRing Mulike','MultiRing MultiGeV Elike Other'])

# First plot: reconstructed momentum distributions
fig, axes = plt.subplots(nrows=4, ncols=4)
fig.tight_layout(h_pad=1)
axis = axes.flat

for it in range(14):
	cc = np.extract((abs(mode)<30) & (itype==it), evis)
	nc = np.extract((abs(mode)>30) & (itype==it), evis)
	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
	weights = [np.ones(np.size(cc[nu==12]))*norm, np.ones(np.size(cc[nu==-12]))*norm, np.ones(np.size(cc[abs(nu)==14]))*norm, np.ones(np.size(cc[abs(nu)==16]))*norm, np.ones(np.size(nc))*norm]
	if it<2:
		bins = sge_ebins
	elif it==2:
		bins = sgsrpi0ebins
	elif it<6:
		bins = sgm_ebins
	elif it==6:
		bins = sgmrpi0ebins
	elif it<9:
		bins = mge_ebins
	elif it==9:
		bins = mgm_ebins
	elif it<12:
		bins = mre0_ebins
	elif it==12:
		bins = mrm_ebins
	elif it==13:
		bins = mre1_ebins
	axis[it].hist(series, bins, weights=weights, stacked=True)
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Momentum (GeV)")
	axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})


plt.show()
# plt.savefig('RecoMomentum_Hists_Unosc.png')
plt.clf()

# First plot: reconstructed zenith distributions
fig, axes = plt.subplots(nrows=4, ncols=4)
fig.tight_layout(h_pad=1)
axis = axes.flat

for it in range(14):
	cc = np.extract((abs(mode)<30) & (itype==it), cz)
	nc = np.extract((abs(mode)>30) & (itype==it), cz)
	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
	weights = [np.ones(np.size(cc[nu==12]))*norm, np.ones(np.size(cc[nu==-12]))*norm, np.ones(np.size(cc[abs(nu)==14]))*norm, np.ones(np.size(cc[abs(nu)==16]))*norm, np.ones(np.size(nc))*norm]
	if it==2 or it==1 or it==6 or it==5:
		bins = z1bins
	else:
		bins = z10bins
	axis[it].hist(series, bins, weights=weights, stacked=True)
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Zenith (GeV)")
	axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})


plt.show()
# plt.savefig('RecoZenith_Hists_Unosc.png')
plt.clf()

