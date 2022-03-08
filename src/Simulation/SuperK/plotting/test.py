import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import applications as ap
from math import pi

# with h5py.File('../../../../../MC/SuperK/data/output/SK/combined.hdf5', 'r') as hf:
# with h5py.File('../data/output/combined_at50deg.hdf5', 'r') as hf:
# with h5py.File('../data/output/combined_last+h.hdf5', 'r') as hf:
with h5py.File('../data/output/combined_last.hdf5', 'r') as hf:
	evis = np.array(hf['evis'][()])
	cx = np.array(hf['recodirX'][()])
	cy = np.array(hf['recodirY'][()])
	cz = np.array(hf['recodirZ'][()])
	# cz = np.array(hf['dirnuZ'][()])
	dz = np.array(hf['dirnuZ'][()])
	dx = np.array(hf['dirnuX'][()])
	dy = np.array(hf['dirnuY'][()])
	mode = np.array(hf['mode'][()])
	ipnu = np.array(hf['ipnu'][()])
	nring = np.array(hf['nring'][()])
	pnu = np.array(hf['pnu'][()])
	oscw = np.array(hf['weightOsc_SKbest'][()])
	weight = np.array(hf['weightReco'][()])
	itype = np.array(hf['itype'][()])
'''
for i, (x, y, z, it) in enumerate(zip(dx, dy, dz, itype)):
	if it == 15:
		ang = np.random.normal(0, 2, 1)
		ang = ang*pi/180.
		u = ap.RndVector()
		di = ap.RodRot(np.array([x,y,z]),u,ang)
		cz[i] = di[2]
'''
dot = dx*cx+dy*cy+dz*cz

wght = oscw*weight

zbins = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
zbins_up = np.array([-0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05])
sge_bins = np.array([0.175, 0.325, 0.515, 0.815, 1.115])
srpi_bins = sge_bins
sgm_bins = sge_bins
mrpi_bins = np.array([0.125, 0.2, 0.325, 0.515, 0.98])

# Exposure normalizations
entries = 39891.6
cond = (itype>-1)*(itype<16)
total = np.sum(wght[cond])
norm = entries / total

print('Number of events from SK paper ', entries)
print('SK Simulation of ', 5326/norm,' days')
print('SK Simulation of ', 5326/norm/365.25,' years')

cond = (itype>-1)*(itype<16)
total = np.sum(wght[cond])
norm = entries / total

# Energy bining
sge_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
sgm_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
sgsrpi0ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
sgmrpi0ebins = np.array([0.1, 0.15, 0.25, 0.4, 0.63, 1.33])
mge_ebins = np.array([1.3, 2.5, 5., 10., 100.])
mgm_ebins = np.array([1.3, 3.0, 100.])
mre_ebins = np.array([1.3, 2.5, 5.0, 100.])
mrm_ebins = np.array([0.6, 1.3, 2.5, 5., 100.])
mro_ebins = np.array([1.3, 2.5, 5.0, 10., 100.])
pcs_ebins = np.array([1.0, 10., 1.0e5])
pct_ebins = np.array([1.0, 10.0, 50., 1.0e3, 1.0e5])
ums_ebins = np.array([1.6, 10.0, 50., 1.0e5])
um_ebins  = np.array([1.6, 1.0e5])

# Zenith bining
z10bins = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
z10bins_up = np.array([-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
z1bins = np.array([-1, 1.0])

# Names of SK FC samples
SampleLabel = np.array(['SingleRingSubGeV Elike 0dcy-e','SingleRing SubGeV Elike 1dcy-e','SingleRing Pi0like','SingleRing SubGeV Mulike 0dcy-e',
	'SingleRing SubGeV Mulike 1dcy-e','SingleRing SubGeV Mulike 2dcy-e','Double Ring Pi0like','SingleRing MultiGeV Elike Nulike','SingleRing MultiGeV Elike NuBarlike',
	'SingleRing MultiGeV Mulike','MultiRing MultiGeV Elike Nulike','MultiRing MultiGeV Elike NuBarlike','MultiRing Mulike','MultiRing MultiGeV Elike Other','PC Stopping',
	'PC Through-going','Upward-going Muons Stopping','Upward-going Muons non-Showering','Upward-going Muons Showering'])

# Second plot: reconstructed zenith distributions
fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(12,12))
fig.tight_layout(h_pad=5)
fig.subplots_adjust(bottom=0.1, top=0.9)
axis = axes.flat

pull = (pnu-evis) / pnu
col = ['green','red','blue','magenta','black']

for i in range(16):
	series = np.extract((itype==i) & (np.abs(pull)<1), pull)
	# series = np.extract((itype==i), dot)
	# weights = np.extract((itype==i), wght)
	# nrng = np.extract((itype==i), nring)
	# weights = np.extract((itype==i) & (np.abs(pull)<2), wght)
	# weights = [np.extract((itype==i) & (abs(mode)<30) & (ipnu==12), weight), np.extract((itype==i) & (abs(mode)<30) & (ipnu==-12), weight), np.extract((itype==i) & (abs(mode)<30) & (abs(ipnu)==14), weight), np.extract((itype==i) & (abs(mode)<30) & (abs(ipnu)==16), weight), np.extract((itype==i) & (abs(mode)>30), weight)]
	# enu = [np.extract((itype==i) & (abs(mode)<30) & (ipnu==12), pnu), np.extract((itype==i) & (abs(mode)<30) & (ipnu==-12), pnu), np.extract((itype==i) & (abs(mode)<30) & (abs(ipnu)==14), pnu), np.extract((itype==i) & (abs(mode)<30) & (abs(ipnu)==16), pnu), np.extract((itype==i) & (abs(mode)>30), pnu)]
	# axis[i].set_xlim([-2,2])
	# axis[i].hist(series, bins=100, weights=weights)
	#axis[i].hist(series, bins=100)
	axis[i].hist(series, bins=50)
	'''
	axis[i].semilogx()
	for n in range(5):
		axis[i].scatter(enu[n],weights[n], c=col[n])
	axis[i].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right')
	'''
	axis[i].set_title(SampleLabel[i])

plt.show()

# for it in range(16):
# 	cc = np.extract((abs(mode)<30) & (itype==it), cz)
# 	wcc = np.extract((abs(mode)<30) & (itype==it), wght)
# 	wnc = np.extract((abs(mode)>30) & (itype==it), wght)
# 	nc = np.extract((abs(mode)>30) & (itype==it), cz)
# 	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
# 	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
# 	weights = [wcc[nu==12]*norm, wcc[nu==-12]*norm, wcc[abs(nu)==14]*norm, wcc[abs(nu)==16]*norm, wnc*norm]
# 	if it==2 or it==1 or it==6 or it==5:
# 		bins = z1bins
# 		skflag = 1
# 	else:
# 		skflag = 1
# 		if it<16:
# 			bins = z10bins
# 		else:
# 			bins = z10bins_up
# 	axis[it].hist(series, bins, weights=weights, stacked=True)
# 	axis[it].set_title(SampleLabel[it])
# 	axis[it].set_xlabel("Reco. Zenith (GeV)")
# 	ymin, ymax = axis[it].get_ylim()
# 	axis[it].set_ylim([0,1.5*ymax])
# 	if skflag: axis[it].legend(labels = ('SK official MC', 'CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
# 	else: axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})

# plt.show()
# plt.savefig('../figs/RecoZenith_Hists_Osc.png')
plt.clf()
