import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import applications as ap
from math import pi

with h5py.File('../data/output/SK_Htag/combined.hdf5', 'r') as hf:
#with h5py.File('../del', 'r') as hf:
	evis = np.array(hf['evis'][()])
	cz = np.array(hf['recodirZ'][()])
	# cz = np.array(hf['dirnuZ'][()])
	dz = np.array(hf['dirnuZ'][()])
	dx = np.array(hf['dirnuX'][()])
	dy = np.array(hf['dirnuY'][()])
	mode = np.array(hf['mode'][()])
	ipnu = np.array(hf['ipnu'][()])
	pnu = np.array(hf['pnu'][()])
	oscw = np.array(hf['weightOsc_SKbest'][()])
	weight = np.array(hf['weightReco'][()])
	itype = np.array(hf['itype'][()])

condition1 = (itype<16) * (itype>-1) 
condition2 = (itype<18) * (itype>15) * (evis>1)
condition = (condition2 + condition1) 
evis=evis[condition]
cz=cz[condition]
dz=dz[condition]
mode=mode[condition]
ipnu=ipnu[condition]
pnu=pnu[condition]
weight=weight[condition]
oscw=oscw[condition]
itype=itype[condition]

wght = oscw*weight

# print(sk_zenith)
zbins = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
zbins_up = np.array([-0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05])
sge_bins = np.array([0.175, 0.325, 0.515, 0.815, 1.115])
srpi_bins = sge_bins
sgm_bins = sge_bins
mrpi_bins = np.array([0.125, 0.2, 0.325, 0.515, 0.98])

# Exposure normalizations
# entries = np.sum(sksge0) + np.sum(sksge1) + np.sum(sksrpi) + np.sum(sksgm0) + np.sum(sksgm1) + np.sum(sksgm2) + np.sum(skmrpi) + np.sum(skmge) + np.sum(skmgeb) + np.sum(skmgm) + np.sum(skmre) + np.sum(skmreb) + np.sum(skmrm) + np.sum(skmro) + np.sum(pcstop) + np.sum(pcthru) + np.sum(umstop) + np.sum(umnons) + np.sum(umshow)
norm = 269*365.25
# print('Number of events from SK paper ', entries)
print('SK Simulation of ', norm,' days')
print('SK Simulation of ', norm/365.25,' years')
cond = (itype>-1)*(itype<18)
total = np.sum(wght[cond])
norm = 36437.9607 / total

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
SampleLabel = np.array(['SingleRing SubGeV NuElike','SingleRing SubGeV Elike','SingleRing SubGeV NuEBarlike','SingleRing Pi0like','SingleRing SubGeV NuMulike',
	'SingleRing SubGeV NuMuBarlike','Double Ring Pi0like','SingleRing MultiGeV NuElike','SingleRing MultiGeV Elike', 'SingleRing MultiGeV NuEBarlike',
	'SingleRing MultiGeV NuMulike','SingleRing MultiGeV NuMuBarlike','MultiRing MultiGeV Elike Nulike','MultiRing MultiGeV Elike NuBarlike','MultiRing Mulike',
	'MultiRing MultiGeV Elike Other','PC Stopping','PC Through-going','Upward-going Muons Stopping','Upward-going Muons non-Showering','Upward-going Muons Showering'])

# Second plot: reconstructed zenith distributions
fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(16,12))
fig.tight_layout(h_pad=6)
fig.subplots_adjust(bottom=0.1, top=0.9)
axis = axes.flat

# for it in range(19):
for it in range(18):
	cc = np.extract((abs(mode)<30) & (itype==it), cz)
	print(np.sum(cc))
	wcc = np.extract((abs(mode)<30) & (itype==it), wght)
	wnc = np.extract((abs(mode)>30) & (itype==it), wght)
	nc = np.extract((abs(mode)>30) & (itype==it), cz)
	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
	weights = [wcc[nu==12]*norm, wcc[nu==-12]*norm, wcc[abs(nu)==14]*norm, wcc[abs(nu)==16]*norm, wnc*norm]
	if it==3 or it==0 or it==6:
		bins = z1bins
		skflag = 0
	else:
		skflag = 0
		if it<18:
			bins = z10bins
		else:
			bins = z10bins_up
	axis[it].hist(series, bins, weights=weights, stacked=True)
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Zenith (GeV)")
	ymin, ymax = axis[it].get_ylim()
	axis[it].set_ylim([0,1.5*ymax])
	if skflag: axis[it].legend(labels = ('SK official MC', 'CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
	else: axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})

plt.show()
# plt.savefig('../figs/RecoZenith_Hists_Osc.png')
# plt.clf()

# Third plot: Ratio plots of the previous
fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(12,12))
fig.tight_layout(h_pad=5)
fig.subplots_adjust(bottom=0.1, top=0.9)
axis = axes.flat

for it in range(16):
	if it==3 or it==0 or it==6:
		bins = z1bins
		bins_mp = np.zeros(1)
		nbins=1
	else:
		bins = z10bins
		bins_mp = zbins
		nbins=10
	sim_hist , dummybins = np.histogram(cz[itype==it], bins=bins, weights=wght[itype==it]*norm)
	sim_entries , dummybins = np.histogram(cz[itype==it], bins=bins)
	# print(sim_entries)
	y_error = np.divide(sim_hist,np.sqrt(sim_entries))
	axis[it].plot(bins,np.zeros(nbins+1))
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Zenith (GeV)")
	axis[it].set_ylabel("SKSim / SKOfficial")
	axis[it].errorbar(bins_mp, np.ones(nbins)-np.divide(sim_hist,sk_zenith[it]), yerr=np.divide(y_error,sk_zenith[it]), fmt='o')
	ymin, ymax = axis[it].get_ylim()
	axis[it].set_ylim([-0.25,0.25])

plt.show()
# plt.savefig('../figs/RecoZenith_Ratios_Osc.png')
# plt.clf()



fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
fig.tight_layout(h_pad=3)
fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9)
axis = axes.flat
bins = np.geomspace(0.1, 1e5, num=51)

axis[0].hist((pnu[itype < 2 & (abs(ipnu)!=16)]) ,bins,  weights=weight[itype < 2 & (abs(ipnu)!=16)], ec='red', fc='none', histtype='step')
axis[0].hist((pnu[(itype > 6) & (itype < 9) & (abs(ipnu)!=16)]), bins, weights=weight[(itype > 6) & (itype < 9) & (abs(ipnu)!=16)], ec='blue', fc='none', histtype='step')
axis[0].hist((pnu[(itype>9) & (itype!=12)]),bins, weights=weight[(itype>9) & (itype!=12)], ec='green', fc='none' , histtype='step')
axis[0].set_ylabel('Events (a.u.)')
axis[0].set_xlabel('Neutrino Energy (GeV)')
axis[0].semilogx()
axis[0].legend(labels = ('SuGeV e-like', 'MultiGeV e-like', 'MultiRing e-like'), loc = 'upper right')

axis[1].hist((pnu[(itype>2) & (itype<6) & (abs(ipnu)!=16)]) ,bins, weights=weight[(itype>2) & (itype<6) & (abs(ipnu)!=16)], ec='red', fc='none', histtype='step')
axis[1].hist((pnu[(itype == 9) & (abs(ipnu)!=16)]), bins, weights=weight[(itype == 9) & (abs(ipnu)!=16)], ec='blue', fc='none', histtype='step')
axis[1].hist((pnu[(itype == 12) & (abs(ipnu)!=16)]),bins, weights=weight[(itype == 12) & (abs(ipnu)!=16)], ec='green', fc='none' , histtype='step')
axis[1].set_ylabel('Events (a.u.)')
axis[1].set_xlabel('Neutrino Energy (GeV)')
axis[1].semilogx()
axis[1].legend(labels = ('SuGeV mu-like', 'MultiGeV mu-like', 'MultiRing mu-like'), loc = 'upper right')

axis[2].hist((pnu[(itype==14) & (abs(ipnu)!=16)]) ,bins, weights=weight[(itype==14) & (abs(ipnu)!=16)], ec='red', fc='none', histtype='step')
axis[2].hist((pnu[(itype==15) & (abs(ipnu)!=16)]) ,bins, weights=weight[(itype==15) & (abs(ipnu)!=16)], ec='blue', fc='none', histtype='step')
axis[2].set_ylabel('Events (a.u.)')
axis[2].set_xlabel('Neutrino Energy (GeV)')
axis[2].semilogx()
axis[2].legend(labels = ('PC Stop', 'PC Thru'), loc = 'upper right')

# plt.show()
plt.savefig('../figs/NeutrinoEnergy_Samples_Unosc.png')
plt.clf()


# First plot: reconstructed momentum distributions
fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(12,12))
fig.tight_layout(h_pad=5)
fig.subplots_adjust(bottom=0.1, top=0.9)
axis = axes.flat

axis[1].plot(sge_bins, sksge1, 'ko')
axis[2].plot(srpi_bins, sksrpi, 'ko')
axis[5].plot(sgm_bins, sksgm2, 'ko')
axis[6].plot(mrpi_bins, skmrpi, 'ko')
skflag = 0

# for it in range(19):
for it in range(16):
	cc = np.extract((abs(mode)<30) & (itype==it), evis)
	wcc = np.extract((abs(mode)<30) & (itype==it), wght)
	nc = np.extract((abs(mode)>30) & (itype==it), evis)
	wnc = np.extract((abs(mode)>30) & (itype==it), wght)
	nu = np.extract((abs(mode)<30) & (itype==it), ipnu)
	series = [cc[nu==12], cc[nu==-12], cc[abs(nu)==14], cc[abs(nu)==16], nc]
	weights = [wcc[nu==12]*norm, wcc[nu==-12]*norm, wcc[abs(nu)==14]*norm, wcc[abs(nu)==16]*norm, wnc*norm]
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
		bins = mre_ebins
	elif it==12:
		bins = mrm_ebins
	elif it==13:
		bins = mro_ebins
	elif it==14:
		bins = pcs_ebins
	elif it==15:
		bins = pct_ebins
	# elif it==16:
	# 	bins = ums_ebins
	# elif it>16:
	# 	bins = um_ebins
	axis[it].hist(series, bins, weights=weights, stacked=True)
	axis[it].set_title(SampleLabel[it])
	axis[it].set_xlabel("Reco. Momentum (GeV)")
	if skflag: axis[it].legend(labels = ('SK official MC', 'CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
	else: axis[it].legend(labels = ('CCnue', 'CCnuebar', 'CCnumu', 'CCnutau', 'NC'), loc = 'upper right', ncol=2, prop={'size': 6})
# plt.show()
plt.savefig('../figs/RecoMomentum_Hists_Osc.png')
