import numpy as np
from plotting import cornerPlot, cornerPlotBothO
import sys

experiment = str(sys.argv[1])
filename = str(sys.argv[2])
if experiment=='SuperK':
	NT23 = 25
	SqT23_min = 0.305
	SqT23_max = 0.705
	SqT23 = np.linspace(SqT23_min, SqT23_max, NT23, endpoint = True)
	SqT23BF = SqT23[16]
	NDM31 = 25
	DM31_max = 3.0e-3
	DM31_min = 2.0e-3
	DM31 = np.linspace(DM31_min, DM31_max, NDM31, endpoint = True)
	DM31BF = DM31[12]

	data = 'SK/X2_stat_25_25_NO_NO.txt'
	linesZ = open(data, "r").readlines()
	xi = np.zeros(NT23*NDM31)
	for i, line in enumerate(linesZ):
		xi[i] = float(line.split()[0])

	data = 'SK/X2_stat_25_25_NO_IO.txt'
	linesZ = open(data, "r").readlines()
	xii = np.zeros(NT23*NDM31)
	for i, line in enumerate(linesZ):
		xii[i] = float(line.split()[0])

	# cornerPlot(SqT23, DM31, xi)
	# cornerPlot(SqT23, DM31, xii)
	cornerPlotBothO(SqT23, DM31, xi, xii)


elif experiment=='IceCube':
	NT23=25
	NDM31=20
	SqT23 = np.zeros(NT23)
	DM31 = np.zeros(NDM31)
	data = 'IC/Xi_ST23_DM31_NoSys.dat'
	data = 'SK/X2_stat_25_20_NO_NO.txt'
	linesZ = open(data, "r").readlines()
	xi = np.zeros(NT23*NDM31)
	for i, line in enumerate(linesZ):
		x, y, z = line.split( )
		xi[i] = float(z)
		if i%NDM31==0:
			SqT23[i//NDM31] = float(x)
		if i<NDM31:
			DM31[i] = float(y)
	data = 'SK/X2_stat_25_20_NO_IO.txt'
	linesZ = open(data, "r").readlines()
	xii = np.zeros(NT23*NDM31)
	for i, line in enumerate(linesZ):
		x, y, z = line.split( )
		xii[i] = float(z)
	cornerPlotBothO(SqT23, DM31, xi, xii)
	# cornerPlot(SqT23,DM31,xi)
