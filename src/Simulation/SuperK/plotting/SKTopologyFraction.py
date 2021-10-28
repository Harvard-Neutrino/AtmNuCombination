import matplotlib.pyplot as plt
import numpy as np

# Initialize
loge = np.zeros(60)
fce = np.zeros(60)
fcm = np.zeros(60)
pcs = np.zeros(60)
pct = np.zeros(60)
ums = np.zeros(60)
umt = np.zeros(60)
umsh = np.zeros(60)

# Acquiring digitized data
with open('../lib/SKTopologyFraction.dat') as f:
    lines = f.readlines()
    for i,l in enumerate(lines):
        loge[i], dummy, fce[i], fcm[i], pcs[i], pct[i], ums[i], umt[i], umsh[i] = l.split( )

# Plotting samples with no flavour information
plt.plot([],[],color='m', label=r"FC $\nu_{e}+\bar{\nu_e}$", linewidth=5)
plt.plot([],[],color='c', label=r"FC $\nu_{\mu}+\bar{\nu_\mu}$", linewidth=5)
plt.plot([],[],color='r', label="PC Stop", linewidth=5)
plt.plot([],[],color='b', label='PC Through', linewidth=5)
plt.plot([],[],color='g', label=r"Up$\mu$ Stop", linewidth=5)
plt.plot([],[],color='y', label=r"Up$\mu$ Through", linewidth=5)
plt.plot([],[],color='coral', label=r"Up$\mu$ Showering", linewidth=5)
plt.stackplot(loge, fce,fcm,pcs,pct,ums,umt,umsh, colors=['m','c','r','b','g','y','coral'])
plt.xlabel(r"$\log_{10}$($E_\nu$ [GeV])")
plt.ylabel('SK sample fraction')
plt.title('Extracted from unoscillated flux -- No flavour information')
plt.legend()
plt.show()
plt.savefig('../figs/FC_PC_UpMu.png')

# Plotting samples with flavour information

# Add a bit of flavour
## electronic
nue = fce + 0.116*pcs + 0.009*pct + 0.011*ums + 0.003*umt + 0.001*umsh
fc_nue   = fce / nue
pcs_nue  = 0.116*pcs / nue
pct_nue  = 0.009*pct / nue
ums_nue  = 0.011*ums / nue
umt_nue  = 0.003*umt / nue
umsh_nue = 0.001*umsh / nue
plt.plot([],[],color='m', label=r"FC", linewidth=5)
plt.plot([],[],color='r', label=r"PC Stop", linewidth=5)
plt.plot([],[],color='b', label=r"PC Through", linewidth=5)
plt.plot([],[],color='g', label=r"Up$\mu$ Stop", linewidth=5)
plt.plot([],[],color='y', label=r"Up$\mu$ Through", linewidth=5)
plt.plot([],[],color='coral', label=r"Up$\mu$ Showering", linewidth=5)
plt.stackplot(loge, fc_nue,pcs_nue,pct_nue,ums_nue,umt_nue,umsh_nue, colors=['m','r','b','g','y','coral'])
plt.xlabel(r"$\log_{10}$($E_\nu$ [GeV])")
plt.ylabel('SK sample fraction')
plt.title(r"Sample by sample distribution of $\nu_{e}+\bar{\nu_e}$")
plt.legend(loc='lower left')
plt.show()
plt.savefig('../figs/FC_PC_UpMu-nue.png')
## muonic
numu = fcm + 0.829*pcs + 0.978*pct + 0.986*ums + 0.996*umt + 0.998*umsh
fc_numu   = fcm / numu
pcs_numu  = 0.829*pcs / numu
pct_numu  = 0.978*pct / numu
ums_numu  = 0.986*ums / numu
umt_numu  = 0.996*umt / numu
umsh_numu = 0.998*umsh / numu
plt.plot([],[],color='m', label=r"FC", linewidth=5)
plt.plot([],[],color='r', label=r"PC Stop", linewidth=5)
plt.plot([],[],color='b', label=r"PC Through", linewidth=5)
plt.plot([],[],color='g', label=r"Up$\mu$ Stop", linewidth=5)
plt.plot([],[],color='y', label=r"Up$\mu$ Through", linewidth=5)
plt.plot([],[],color='coral', label=r"Up$\mu$ Showering", linewidth=5)
plt.stackplot(loge, fc_numu,pcs_numu,pct_numu,ums_numu,umt_numu,umsh_numu, colors=['m','r','b','g','y','coral'])
plt.xlabel(r"$\log_{10}$($E_\nu$ [GeV])")
plt.ylabel('SK sample fraction')
plt.title(r"Sample by sample distribution of $\nu_{\mu}+\bar{\nu_\mu}$")
plt.legend(loc='lower left')
plt.show()
plt.savefig('../figs/FC_PC_UpMu-numu.png')
## tauonic
nut = 0.0057*(fce+fcm) + 0.01*pcs + 0.007*pct + 0.0*ums + 0.0*umt + 0.0*umsh
nut[nut==0] = 1.
fc_nut   = 0.0057*(fce+fcm) / nut
pcs_nut  = 0.01*pcs / nut
pct_nut  = 0.007*pct / nut
ums_nut  = 0.0*ums / nut
umt_nut  = 0.0*umt / nut
umsh_nut = 0.0*umsh / nut
plt.plot([],[],color='m', label=r"FC", linewidth=5)
plt.plot([],[],color='r', label=r"PC Stop", linewidth=5)
plt.plot([],[],color='b', label=r"PC Through", linewidth=5)
plt.plot([],[],color='g', label=r"Up$\mu$ Stop", linewidth=5)
plt.plot([],[],color='y', label=r"Up$\mu$ Through", linewidth=5)
plt.plot([],[],color='coral', label=r"Up$\mu$ Showering", linewidth=5)
plt.stackplot(loge, fc_nut,pcs_nut,pct_nut,ums_nut,umt_nut,umsh_nut, colors=['m','r','b','g','y','coral'])
plt.xlabel(r"$\log_{10}$($E_\nu$ [GeV])")
plt.ylabel('SK sample fraction')
plt.title(r"Sample by sample distribution of $\nu_{\tau}+\bar{\nu_\tau}$")
plt.legend(loc='lower left')
plt.show()
plt.savefig('../figs/FC_PC_UpMu-nutau.png')
## NC
nc = 0.12*(fce+fcm) + 0.045*pcs + 0.006*pct + 0.003*ums + 0.001*umt + 0.001*umsh
fc_nc   = 0.12*(fce+fcm) / nc
pcs_nc  = 0.045*pcs / nc
pct_nc  = 0.006*pct / nc
ums_nc  = 0.003*ums / nc
umt_nc  = 0.001*umt / nc
umsh_nc = 0.001*umsh / nc
plt.plot([],[],color='m', label=r"FC", linewidth=5)
plt.plot([],[],color='r', label=r"PC Stop", linewidth=5)
plt.plot([],[],color='b', label=r"PC Through", linewidth=5)
plt.plot([],[],color='g', label=r"Up$\mu$ Stop", linewidth=5)
plt.plot([],[],color='y', label=r"Up$\mu$ Through", linewidth=5)
plt.plot([],[],color='coral', label=r"Up$\mu$ Showering", linewidth=5)
plt.stackplot(loge, fc_nc,pcs_nc,pct_nc,ums_nc,umt_nc,umsh_nc, colors=['m','r','b','g','y','coral'])
plt.xlabel(r"$\log_{10}$($E_\nu$ [GeV])")
plt.ylabel('SK sample fraction')
plt.title(r"Sample by sample distribution of NC events")
plt.legend(loc='lower left')
plt.show()
plt.savefig('../figs/FC_PC_UpMu-nc.png')
