import numpy as np
import pandas as pd
from plotting import cornerPlot, cornerPlotBothO
import sys


def X2reader(filename, filename2='None'):

	df = pd.read_csv(filename, sep=' ')
	rawOscPar = {}
	oscPar = {}
	for i in range(df.columns.size-1):
		rawOscPar[i] = df.values[:,i]
		oscPar[i] = np.unique(rawOscPar[i])
	chi2 = np.array(df["X2"])

	if df.columns.size == 3:
		cornerPlot(oscPar[0], oscPar[1], chi2)
	elif df.columns.size == 4:
		X2 = np.array([])
		for q in oscPar[0]:
			for m in oscPar[1]:
				cut = (rawOscPar[0]==q) * (rawOscPar[1]==m)
				chi2_slice = chi2[cut]
				marg = np.amin(chi2_slice)
				X2 = np.append(X2,marg)
		cornerPlot(oscPar[0],oscPar[1],X2)
		X2 = np.array([])
		for q in oscPar[0]:
			for m in oscPar[2]:
				cut = (rawOscPar[0]==q) * (rawOscPar[2]==m)
				chi2_slice = chi2[cut]
				marg = np.amin(chi2_slice)
				X2 = np.append(X2,marg)
		cornerPlot(oscPar[0],oscPar[2],X2)

	if filename2!='None':
		df2 = pd.read_csv(filename2, sep=' ')
		chi2I = np.array(df2["X2"])

		if df.columns.size == 3:
			cornerPlotBothO(oscPar[0],oscPar[1], chi2, chi2I)
		elif df.columns.size == 4:
			X2I = np.array([])
			for q in oscPar[0]:
				for m in oscPar[1]:
					cut = (rawOscPar[0]==q) * (rawOscPar[1]==m)
					chi2_slice = chi2I[cut]
					marg = np.amin(chi2_slice)
					X2I = np.append(X2,marg)
			cornerPlotBothO(oscPar[0],oscPar[1],X2,X2I)
			X2I = np.array([])
			for q in oscPar[0]:
				for m in oscPar[2]:
					cut = (rawOscPar[0]==q) * (rawOscPar[2]==m)
					chi2_slice = chi2I[cut]
					marg = np.amin(chi2_slice)
					X2I = np.append(X2I,marg)
			cornerPlotBothO(oscPar[0],oscPar[2],X2,X2I)

def X2readerBoth(filename, filename2):

	df = pd.read_csv(filename, sep=' ')
	rawOscPar = {}
	oscPar = {}
	for i in range(df.columns.size-1):
		rawOscPar[i] = df.values[:,i]
		oscPar[i] = np.unique(rawOscPar[i])
	chi2 = np.array(df["X2"])

	df2 = pd.read_csv(filename2, sep=' ')
	chi2I = np.array(df2["X2"])
	if df.columns.size == 3:
		cornerPlotBothO(oscPar[0],oscPar[1], chi2, chi2I)
	elif df.columns.size == 4:
		X2I = np.array([])
		for q in oscPar[0]:
			for m in oscPar[1]:
				cut = (rawOscPar[0]==q) * (rawOscPar[1]==m)
				chi2_slice = chi2I[cut]
				marg = np.amin(chi2_slice)
				X2I = np.append(X2,marg)
		cornerPlotBothO(oscPar[0],oscPar[1],X2,X2I)
		X2I = np.array([])
		for q in oscPar[0]:
			for m in oscPar[2]:
				cut = (rawOscPar[0]==q) * (rawOscPar[2]==m)
				chi2_slice = chi2I[cut]
				marg = np.amin(chi2_slice)
				X2I = np.append(X2I,marg)
		cornerPlotBothO(oscPar[0],oscPar[2],X2,X2I)

	# Combination...
	chi2 = chi2+chi2I
	if df.columns.size == 3:
		cornerPlotBothO(oscPar[0],oscPar[1], chi2, title='SuperK and DeepCore combination')
	elif df.columns.size == 4:
		X2 = np.array([])
		for q in oscPar[0]:
			for m in oscPar[1]:
				cut = (rawOscPar[0]==q) * (rawOscPar[1]==m)
				chi2_slice = chi2[cut]
				marg = np.amin(chi2_slice)
				X2 = np.append(X2,marg)
		cornerPlot(oscPar[0],oscPar[1],X2, title='SuperK and DeepCore combination')
		X2 = np.array([])
		for q in oscPar[0]:
			for m in oscPar[2]:
				cut = (rawOscPar[0]==q) * (rawOscPar[2]==m)
				chi2_slice = chi2[cut]
				marg = np.amin(chi2_slice)
				X2 = np.append(X2,marg)
		cornerPlot(oscPar[0],oscPar[2],X2, title='SuperK and DeepCore combination')




experiment = str(sys.argv[1])
filename = str(sys.argv[2])
if len(sys.argv)==3 and (experiment=='SK' or experiment=='IC'):
	X2reader(filename)
elif len(sys.argv)==4 and (experiment=='SK' or experiment=='IC'):
	filename2 = sys.argv[3]
	X2reader(filename, filename2)
elif experiment=='SK+IC' or experiment='IC+SK':
	if len(sys.argv)<4:
		print('You are missing files')
	elif len(sys.argv)==4:
		print('You are combining IC and SK X2 surfaces, make sure they share common points...')
		filename2 = sys.argv[3]
		X2readerBoth(filename,filename2)