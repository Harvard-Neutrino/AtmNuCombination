import numpy as np
from scipy import interpolate
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager as font_manager
import pandas as pd
import scipy as scp

matplotlib.rcParams.update({'font.size': 21})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})


font1 = font_manager.FontProperties(
                           weight='bold',
                           style='normal', size=16)
font2 = font_manager.FontProperties(
                           style='normal', size=16)

ICsep = pd.read_csv("./../data/ICUp_t23_dm31.dat", delimiter = ' ')
SKsep = pd.read_csv("./../data/SK2020+SKGd_t23_dm31.dat", delimiter = ' ')

theta23 = ICsep['Sin2Theta23']
dm31 = ICsep['Dm231']
chisq1 = ICsep['X2']

theta232 = SKsep["Sin2Theta23"]
dm312 = SKsep["Dm231"]
chisq2 = SKsep["X2"]

X, Y = np.meshgrid(theta23, dm31)
Z1 = chisq1
Z2 = chisq2

t23min = 0.3
t23max = 0.7
m31min = 0.00235
m31max = 0.00265
ogt23step = 0.001 * np.pi
ogm31step = 0.01e-3

Z1interp = scp.interpolate.interp2d(theta23, dm31, Z1, kind = 'cubic')
Z2interp = scp.interpolate.interp2d(theta232, dm312, Z2, kind = 'cubic')
newx = np.arange(t23min, t23max + ogt23step, ogt23step)
newy = np.arange(m31min, m31max + ogm31step, ogm31step)
newz1 = Z1interp(newx, newy)
newz2 = Z2interp(newx, newy)
XX, YY = np.meshgrid(newx, newy)

fig, ax  = plt.subplots(figsize=(13,10))
fig.suptitle("Separate and Combined Sensitivities")
ax.set_xlabel(r"$\sin^2{\theta_{23}}$")
ax.set_ylabel(r"$\Delta m^2_{31}$ [eV^2]")
colors = ['orange', 'cyan', 'darkgreen']
axim = ax.contour(XX,YY,newz1,levels=[4.605, 9.21],colors = ["orchid", "mediumvioletred"], linestyles = ['--', '-.'], \
								linewidths = 2.5, label = "with systematics", alpha = 0.6)
axim2 = ax.contour(XX,YY,newz2, levels=[4.605, 9.21], \
									colors = ["lightsalmon", "peru"], linestyles = ['--', '-.'], \
									linewidths = 2.5, label = "without systematics", alpha = 0.6)
h1,_ = axim.legend_elements()
h2,_ = axim2.legend_elements()
legend1 = ax.legend([h2[0], h2[1], h1[0], h1[1]], [r'SK+SKGd(5 yrs) @ 90% CL',\
  r'SK+SKGd(5 yrs) @ 99% CL', r'IC Upgrade(5 yrs) @ 90% CL',\
  r'IC Upgrade(5 yrs) @ 99% CL'], loc = 2, prop = font2)


sys = pd.read_csv("./../data/XiDM31_T23_Sys.dat", sep = " ")
combtheta23 = sys["Theta23"]
combdm31 = sys["DM31"]
combchisq = sys["ChiSq"]

combX, combY = np.meshgrid(combtheta23, combdm31)
combZ = combchisq


combZinterp = scp.interpolate.interp2d(combtheta23, combdm31, combZ, kind = 'cubic')
combnewx = np.arange(t23min, t23max + ogt23step, ogt23step)
combnewy = np.arange(m31min, m31max + ogm31step, ogm31step)
combnewz = combZinterp(combnewx, combnewy)
combXX, combYY = np.meshgrid(combnewx, combnewy)
colors_comb = ['sienna', 'lawngreen', 'olivedrab']

axim_comb = ax.contour(combXX,combYY,combnewz,levels=[2.71, 3.84, 6.63],colors = colors_comb, \
							linewidths = 3)
h_comb,_ = axim_comb.legend_elements()
legend_comb = ax.legend([h_comb[0], h_comb[1], h_comb[2]], [r'SK+SKGd(5yrs)+ICUp(5yrs) @ 90% CL',\
r'SK+SKGd(5yrs)+ICUp(5yrs) @ 95% CL', r'SK+SKGd(5yrs)+ICUp(5yrs) @ 99% CL'], loc = 4, prop = font2)
plt.gca().add_artist(legend_comb)

plt.gca().add_artist(legend1)
ax.ticklabel_format(style = 'sci', scilimits = (-2,1))
# plt.show()
fig.savefig("./../plots/Combined_vs_Separate", bbox_inches="tight")
plt.close()