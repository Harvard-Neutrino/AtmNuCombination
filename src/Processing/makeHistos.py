import h5py
import numpy as np
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

with h5py.File('../Simulation/SuperK/data/testfcmc.hdf5', 'r') as hf:
	evis = np.array(hf['evis'][()])
	cz = np.array(hf['recodirZ'][()])
	# cz = np.array(hf['dirnuZ'][()])
	mode = np.array(hf['mode'][()])
	ipnu = np.array(hf['ipnu'][()])
	pnu = np.array(hf['pnu'][()])
	oscw = np.array(hf['oscw'][()])
	itype = np.array(hf['itype'][()])

# SK results ugly part
sksge0 = np.array([1035.6394129979037,1060.7966457023063,1054.5073375262057,1054.5073375262057,1048.218029350105,1048.218029350105,1035.6394129979037,1020.9643605870023,985.3249475890988,970.6498951781972])
sksge1 = np.array([287.68768768768774,233.03303303303306,252.8528528528529,242.64264264264273,140.54054054054058])
sksrpi = np.array([236.24454148471617,157.2052401746725,99.56331877729258,55.021834061135394,22.707423580786042])
sksgm0 = np.array([245.96328610776658,253.91081569933965,256.09561665374696,263.3955127465047,272.6350106930998,283.1400878538481,296.86024201041266,308.01295266997624,319.17555925736553,327.74763184949387])
sksgm1 = np.array([668.5995853739145,681.2997626415886,690.4786227203074,724.3277348796681,772.3101883844603,806.1773277649249,843.562779797494,888.0028843553763,934.2066520445871,962.7767929573653])
sksgm2 = np.array([61.17956966097185,150.80060620643945,185.24503657923952,188.68470932907803,96.73129772925876])
skmrpi = np.array([550.7141526868473,585.9461627723738,422.59358480294543,154.6795982980674,28.851195507356465])
skmge  = np.array([52.704929803806266,59.86336873345964,64.4919409435627,78.46347951556764,90.87517865324519,89.08556892083183,73.8708089145405,60.79715248952017,54.33794013663129,44.765545396729905])
skmgeb = np.array([175.59607293127635,192.42636746143063,220.47685834502113,261.43057503506316,290.04207573632544,289.4810659186536,255.82047685834507,206.45161290322585,172.2300140252455,146.98457223001407])
skmgm  = np.array([164.18097683941443,176.7056818585275,183.16508278400414,188.28479117601813,223.04342165994063,286.05867663940734,330.91115660911044,340.74937571748114,331.04488916766275,332.7905132726606])
skmre  = np.array([64.97671416055701,74.22372850094524,86.15760593904182,103.73680084843454,133.95213722506568,133.52561442338725,103.25402314750771,81.8550283580025,65.83252639830312,55.18651726840966])
skmreb = np.array([56.618141744503305,62.18110375968402,74.46064693475759,93.45223108704356,114.8435081578828,110.8147151448112,90.47758594359442,67.2652828740182,56.76136387087723,47.456465839252786])
skmrm  = np.array([144.3181818181818,153.97727272727275,169.3181818181818,189.20454545454544,231.25000000000003,276.1363636363637,294.31818181818187,286.3636363636364,273.8636363636364,264.2045454545455])
skmro  = np.array([136.48648648648648,140.99099099099098,158.55855855855853,183.33333333333331,225.2252252252252,237.38738738738735,207.2072072072072,180.18018018018017,158.55855855855853,146.8468468468468])
zbins = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
sge_bins = np.array([0.175, 0.325, 0.515, 0.815, 1.115])
srpi_bins = sge_bins
sgm_bins = sge_bins
mrpi_bins = np.array([0.125, 0.2, 0.325, 0.515, 0.98])

# Livetime normalizations
# years=14.58
# entries = 10266.1 + 1150.7 + 2824.3 + 8008.7 + 687.0 + 571.8 + 1728.4 + 671.3 + 2193.7 + 2573.8 + 915.5 + 773.8 + 2294.0 + 1772.6
entries = np.sum(sksge0) + np.sum(sksge1) + np.sum(sksrpi) + np.sum(sksgm0) + np.sum(sksgm1) + np.sum(sksgm2) + np.sum(skmrpi) + np.sum(skmge) + np.sum(skmgeb) + np.sum(skmgm) + np.sum(skmre) + np.sum(skmreb) + np.sum(skmrm) + np.sum(skmro)
# livetime = years*365.25
total = np.sum(oscw)
# norm = livetime / (total / 8.3)
norm = entries / total

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

axis[1].plot(sge_bins, sksge1, 'ko')
axis[2].plot(srpi_bins, sksrpi, 'ko')
axis[5].plot(sgm_bins, sksgm2, 'ko')
axis[6].plot(mrpi_bins, skmrpi, 'ko')
skflag = 0

for it in range(14):
	cc = np.extract((abs(mode)<30) & (itype==it), evis)
	wcc = np.extract((abs(mode)<30) & (itype==it), oscw)
	nc = np.extract((abs(mode)>30) & (itype==it), evis)
	wnc = np.extract((abs(mode)>30) & (itype==it), oscw)
	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
	weights = [wcc[nu==12]*norm, wcc[nu==-12]*norm, wcc[abs(nu)==14]*norm, wcc[abs(nu)==16]*norm, wnc]
	# weights = [np.ones(np.size(cc[nu==12]))*norm, np.ones(np.size(cc[nu==-12]))*norm, np.ones(np.size(cc[abs(nu)==14]))*norm, np.ones(np.size(cc[abs(nu)==16]))*norm, np.ones(np.size(nc))*norm]
	if it<2:
		bins = sge_ebins
		if it==1: skflag = 1
	elif it==2:
		bins = sgsrpi0ebins
		skflag = 1
	elif it<6:
		bins = sgm_ebins
		if it==5: skflag = 1
	elif it==6:
		bins = sgmrpi0ebins
		skflag = 1
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
	if skflag: axis[it].legend(labels = ('SK', 'CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
	else: axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})


plt.show()
# plt.savefig('RecoMomentum_Hists_Unosc.png')
plt.clf()

# First plot: reconstructed zenith distributions
fig, axes = plt.subplots(nrows=4, ncols=4)
fig.tight_layout(h_pad=1)
axis = axes.flat

axis[0].plot(zbins, sksge0, 'ko')
axis[3].plot(zbins, sksgm0, 'ko')
axis[4].plot(zbins, sksgm1, 'ko')
axis[7].plot(zbins, skmge, 'ko')
axis[8].plot(zbins, skmgeb, 'ko')
axis[9].plot(zbins, skmgm, 'ko')
axis[10].plot(zbins, skmre, 'ko')
axis[11].plot(zbins, skmreb, 'ko')
axis[12].plot(zbins, skmrm, 'ko')
axis[13].plot(zbins, skmro, 'ko')


for it in range(14):
	cc = np.extract((abs(mode)<30) & (itype==it), cz)
	wcc = np.extract((abs(mode)<30) & (itype==it), oscw)
	wnc = np.extract((abs(mode)>30) & (itype==it), oscw)
	nc = np.extract((abs(mode)>30) & (itype==it), cz)
	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
	weights = [wcc[nu==12]*norm, wcc[nu==-12]*norm, wcc[abs(nu)==14]*norm, wcc[abs(nu)==16]*norm, wnc]
	# weights = [np.ones(np.size(cc[nu==12]))*norm, np.ones(np.size(cc[nu==-12]))*norm, np.ones(np.size(cc[abs(nu)==14]))*norm, np.ones(np.size(cc[abs(nu)==16]))*norm, np.ones(np.size(nc))*norm]
	if it==2 or it==1 or it==6 or it==5:
		bins = z1bins
		skflag = 0
	else:
		skflag = 1
		bins = z10bins
	axis[it].hist(series, bins, weights=weights, stacked=True)
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Zenith (GeV)")
	if skflag: axis[it].legend(labels = ('SK', 'CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
	else: axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})


plt.show()
# plt.savefig('RecoZenith_Hists_Unosc.png')
plt.clf()

fig, axes = plt.subplots(nrows=2, ncols=2)
fig.tight_layout(h_pad=1)
axis = axes.flat


axis[0].hist(np.log10(pnu[itype < 2]) ,50,  weights=oscw[itype < 2], ec='red', fc='none')
axis[0].hist(np.log10(pnu[(itype > 6) & (itype < 9)]), 50, weights=oscw[(itype > 6) & (itype < 9)], ec='blue', fc='none')
axis[0].hist(np.log10(pnu[(itype>9) & (itype!=12)]),50, weights=oscw[(itype>9) & (itype!=12)], ec='green', fc='none' )
axis[0].legend(labels = ('SuGeV e-like', 'MultiGeV e-like', 'MultiRing e-like'), loc = 'upper right')

axis[1].hist(np.log10(pnu[(itype>2) & (itype<6)]) ,50, weights=oscw[(itype>2) & (itype<6)], ec='red', fc='none')
axis[1].hist(np.log10(pnu[(itype == 9)]), 50, weights=oscw[(itype == 9)], ec='blue', fc='none')
axis[1].hist(np.log10(pnu[(itype == 12)]),50, weights=oscw[(itype == 12)], ec='green', fc='none' )
axis[1].legend(labels = ('SuGeV mu-like', 'MultiGeV mu-like', 'MultiRing mu-like'), loc = 'upper right')

plt.show()
