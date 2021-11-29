import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.gridspec import GridSpec
import math
from SimReader import Reader
import sys

experiment=str(sys.argv[1])
infile=str(sys.argv[2])

mc = Reader(experiment,infile)
mc.Binning()

fig, axs = plt.subplots(1, 2, sharey=True, sharex=True, tight_layout=True)

for i in range(2):
	print(mc.EnergyBins[i])
	axs[i].hist(mc.EReco[mc.Sample==i], bins=mc.EnergyBins[i])
	axs[i].semilogx()

plt.show()

c = (mc.EReco>+mc.EnergyBins[i][0]) * (mc.EReco<=mc.EnergyBins[i][1])

binn = mc.EReco[c]
print(binn.size)
