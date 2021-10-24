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
	weight = np.array(hf['weight'][()])
	itype = np.array(hf['itype'][()])

wght = oscw*weight

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
pcstop = np.array([32.0101865840863,37.6299822620335,42.7788228538,48.8706470198338,60.7698438468578,72.1975488867481,73.2669255354263,72.6078567652497,68.6504922889491,59.0422034794323])
pcthru = np.array([142.318630587872,144.836060913942,174.263174236207,237.543676031352,396.310149748476,484.764225959415,407.421390432016,335.286880645815,283.11761953421,274.349615242771])
umstop = np.array([89.3464730165499,94.7406333982608,97.4104898230094,110.435014099743,112.014956207079,126.675315633067,145.148447823725,177.79094797731,226.237205733373,270.324288056356])
umnons = np.array([227.035757057314,279.52369546621,323.99247219846,373.2749358426,422.560136869119,470.242600513259,553.116852010265,651.996578272027,758.87630453379,835.356030795552])
umshow = np.array([69.0617429569478,76.6951866805789,79.2109966882911,94.9460599736878,102.152701537903,123.005761466225,144.711336932359,172.813863811641,211.577371501157,234.562264664519])
sk_zenith = {0:sksge0, 1:np.array([np.sum(sksge1)]), 2:np.array([np.sum(sksrpi)]), 3:sksgm0, 4:sksgm1, 5:np.array([np.sum(sksgm2)]), 6:np.array([np.sum(skmrpi)]), 7:skmge,
8:skmgeb, 9:skmgm, 10:skmre, 11:skmreb, 12:skmrm, 13:skmro, 14:pcstop, 15:pcthru}
# print(sk_zenith)
zbins = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
zbins_up = np.array([-0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05])
sge_bins = np.array([0.175, 0.325, 0.515, 0.815, 1.115])
srpi_bins = sge_bins
sgm_bins = sge_bins
mrpi_bins = np.array([0.125, 0.2, 0.325, 0.515, 0.98])

# Livetime normalizations
entries = np.sum(sksge0) + np.sum(sksge1) + np.sum(sksrpi) + np.sum(sksgm0) + np.sum(sksgm1) + np.sum(sksgm2) + np.sum(skmrpi) + np.sum(skmge) + np.sum(skmgeb) + np.sum(skmgm) + np.sum(skmre) + np.sum(skmreb) + np.sum(skmrm) + np.sum(skmro) + np.sum(pcstop) + np.sum(pcthru) + np.sum(umstop) + np.sum(umnons) + np.sum(umshow)
cond = (itype>-1)*(itype<16)
print(cond.size)
total = np.sum(wght[cond])
norm = entries / total
# years=14.58
# entries = 10266.1 + 1150.7 + 2824.3 + 8008.7 + 687.0 + 571.8 + 1728.4 + 671.3 + 2193.7 + 2573.8 + 915.5 + 773.8 + 2294.0 + 1772.6
# livetime = years*365.25
# norm = livetime / (total / 8.3)

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
pcs_ebins = np.array([1.0, 20., 1.0e5])
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
fig, axes = plt.subplots(nrows=4, ncols=4)
fig.tight_layout(h_pad=1)
axis = axes.flat

axis[0].plot(zbins, sksge0, 'ko')
axis[1].plot(np.zeros(1), [np.sum(sksge1)], 'ko')
axis[2].plot(np.zeros(1), [np.sum(sksrpi)], 'ko')
axis[3].plot(zbins, sksgm0, 'ko')
axis[4].plot(zbins, sksgm1, 'ko')
axis[5].plot(np.zeros(1), np.sum(sksgm2), 'ko')
axis[6].plot(np.zeros(1), np.sum(skmrpi), 'ko')
axis[7].plot(zbins, skmge, 'ko')
axis[8].plot(zbins, skmgeb, 'ko')
axis[9].plot(zbins, skmgm, 'ko')
axis[10].plot(zbins, skmre, 'ko')
axis[11].plot(zbins, skmreb, 'ko')
axis[12].plot(zbins, skmrm, 'ko')
axis[13].plot(zbins, skmro, 'ko')
axis[14].plot(zbins, pcstop, 'ko')
axis[15].plot(zbins, pcthru, 'ko')
# axis[16].plot(zbins_up, umstop, 'ko')
# axis[17].plot(zbins_up, umnons, 'ko')
# axis[18].plot(zbins_up, umshow, 'ko')

# for it in range(19):
for it in range(16):
	cc = np.extract((abs(mode)<30) & (itype==it), cz)
	wcc = np.extract((abs(mode)<30) & (itype==it), wght)
	wnc = np.extract((abs(mode)>30) & (itype==it), wght)
	nc = np.extract((abs(mode)>30) & (itype==it), cz)
	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
	weights = [wcc[nu==12]*norm, wcc[nu==-12]*norm, wcc[abs(nu)==14]*norm, wcc[abs(nu)==16]*norm, wnc]
	# weights = [np.ones(np.size(cc[nu==12]))*norm, np.ones(np.size(cc[nu==-12]))*norm, np.ones(np.size(cc[abs(nu)==14]))*norm, np.ones(np.size(cc[abs(nu)==16]))*norm, np.ones(np.size(nc))*norm]
	if it==2 or it==1 or it==6 or it==5:
		bins = z1bins
		skflag = 1
	else:
		skflag = 1
		if it<16:
			bins = z10bins
		else:
			bins = z10bins_up
	axis[it].hist(series, bins, weights=weights, stacked=True)
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Zenith (GeV)")
	ymin, ymax = axis[it].get_ylim()
	axis[it].set_ylim([0,1.5*ymax])
	if skflag: axis[it].legend(labels = ('SK', 'CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 5})
	else: axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 5})

plt.show()
# plt.savefig('RecoZenith_Hists_Unosc.png')
plt.clf()

# Third plot: Ration plot of the previous
fig, axes = plt.subplots(nrows=4, ncols=5)
fig.tight_layout(h_pad=1)
axis = axes.flat

for it in range(16):
	if it==2 or it==1 or it==6 or it==5:
		bins = z1bins
		bins_mp = np.zeros(1)
		nbins=1
	else:
		bins = z10bins
		bins_mp = zbins
		nbins=10
	sim_hist , dummybins = np.histogram(cz[itype==it], bins=bins, weights=wght[itype==it])
	sim_entries , dummybins = np.histogram(cz[itype==it], bins=bins)
	print(sim_entries)
	y_error = np.divide(sim_hist,np.sqrt(sim_entries))
	axis[it].plot(bins,np.ones(nbins+1))
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Zenith (GeV)")
	axis[it].set_ylabel("SKSim / SKOfficial")
	axis[it].errorbar(bins_mp, np.divide(sim_hist,sk_zenith[it]), yerr=np.divide(y_error,sk_zenith[it]), fmt='o')
	ymin, ymax = axis[it].get_ylim()
	axis[it].set_ylim([0,1.1*ymax])




plt.show()
plt.clf()



# fig, axes = plt.subplots(nrows=2, ncols=2)
# fig.tight_layout(h_pad=1)
# axis = axes.flat

# axis[0].hist(np.log10(pnu[itype < 2]) ,50,  weights=wght[itype < 2], ec='red', fc='none')
# axis[0].hist(np.log10(pnu[(itype > 6) & (itype < 9)]), 50, weights=wght[(itype > 6) & (itype < 9)], ec='blue', fc='none')
# axis[0].hist(np.log10(pnu[(itype>9) & (itype!=12)]),50, weights=wght[(itype>9) & (itype!=12)], ec='green', fc='none' )
# axis[0].legend(labels = ('SuGeV e-like', 'MultiGeV e-like', 'MultiRing e-like'), loc = 'upper right')

# axis[1].hist(np.log10(pnu[(itype>2) & (itype<6)]) ,50, weights=wght[(itype>2) & (itype<6)], ec='red', fc='none')
# axis[1].hist(np.log10(pnu[(itype == 9)]), 50, weights=wght[(itype == 9)], ec='blue', fc='none')
# axis[1].hist(np.log10(pnu[(itype == 12)]),50, weights=wght[(itype == 12)], ec='green', fc='none' )
# axis[1].legend(labels = ('SuGeV mu-like', 'MultiGeV mu-like', 'MultiRing mu-like'), loc = 'upper right')

# axis[2].hist(np.log10(pnu[(itype==14)]) ,50, weights=wght[(itype==14)], ec='red', fc='none')
# axis[2].hist(np.log10(pnu[(itype==15)]) ,50, weights=wght[(itype==15)], ec='blue', fc='none')
# axis[2].legend(labels = ('PC Stop', 'PC Thru'), loc = 'upper right')

# axis[3].hist(np.log10(pnu[(itype==16)]) ,50, weights=wght[(itype==16)], ec='red', fc='none')
# axis[3].hist(np.log10(pnu[(itype==17)]) ,50, weights=wght[(itype==17)], ec='blue', fc='none')
# axis[3].hist(np.log10(pnu[(itype==18)]) ,50, weights=wght[(itype==18)], ec='green', fc='none')
# axis[3].legend(labels = ('UpMu Stop', 'UpMu non-Showering', 'UpMu Showering'), loc = 'upper right')

# plt.show()

# # First plot: reconstructed momentum distributions
# fig, axes = plt.subplots(nrows=4, ncols=5)
# fig.tight_layout(h_pad=1)
# axis = axes.flat

# axis[1].plot(sge_bins, sksge1, 'ko')
# axis[2].plot(srpi_bins, sksrpi, 'ko')
# axis[5].plot(sgm_bins, sksgm2, 'ko')
# axis[6].plot(mrpi_bins, skmrpi, 'ko')
# skflag = 0

# # for it in range(19):
# for it in range(16):
# 	cc = np.extract((abs(mode)<30) & (itype==it), evis)
# 	wcc = np.extract((abs(mode)<30) & (itype==it), wght)
# 	nc = np.extract((abs(mode)>30) & (itype==it), evis)
# 	wnc = np.extract((abs(mode)>30) & (itype==it), wght)
# 	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
# 	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
# 	weights = [wcc[nu==12]*norm, wcc[nu==-12]*norm, wcc[abs(nu)==14]*norm, wcc[abs(nu)==16]*norm, wnc]
# 	# weights = [np.ones(np.size(cc[nu==12]))*norm, np.ones(np.size(cc[nu==-12]))*norm, np.ones(np.size(cc[abs(nu)==14]))*norm, np.ones(np.size(cc[abs(nu)==16]))*norm, np.ones(np.size(nc))*norm]
# 	if it<2:
# 		bins = sge_ebins
# 		if it==1: skflag = 1
# 	elif it==2:
# 		bins = sgsrpi0ebins
# 		skflag = 1
# 	elif it<6:
# 		bins = sgm_ebins
# 		if it==5: skflag = 1
# 	elif it==6:
# 		bins = sgmrpi0ebins
# 		skflag = 1
# 	elif it<9:
# 		bins = mge_ebins
# 	elif it==9:
# 		bins = mgm_ebins
# 	elif it<12:
# 		bins = mre_ebins
# 	elif it==12:
# 		bins = mrm_ebins
# 	elif it==13:
# 		bins = mro_ebins
# 	elif it==14:
# 		bins = pcs_ebins
# 	elif it==15:
# 		bins = pct_ebins
# 	elif it==16:
# 		bins = ums_ebins
# 	elif it>16:
# 		bins = um_ebins
# 	axis[it].hist(series, bins, weights=weights, stacked=True)
# 	axis[it].set_title(SampleLabel[it])
# 	axis[it].set_xlabel("Reco. Momentum (GeV)")
# 	if skflag: axis[it].legend(labels = ('SK', 'CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
# 	else: axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
# plt.show()
# # plt.savefig('RecoMomentum_Hists_Unosc.png')
# plt.clf()
