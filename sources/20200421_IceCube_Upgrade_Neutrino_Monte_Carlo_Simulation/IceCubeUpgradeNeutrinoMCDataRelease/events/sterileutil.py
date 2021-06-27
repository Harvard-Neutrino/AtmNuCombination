import math
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def KronDelta(a, b):
	return (a == b)

def Sigma(lo, hi, func):
	if lo > hi:
		print("Warning")
		print("Sigma: lo > hi encountered")
		return 0
	init = lo
	res = 0
	while init <= hi:
		res += func(init)
		init += 1
	return res

def DoubleSum(lo, hi, func):
	#print("entered DoubleSum function...")
	res = 0
	i = lo
	j = lo
	while j <= hi:
		while i <= j:
			res += func(i,j)
			i += 1
		j += 1
		i = lo
	return res

def TripleCLSum(lo, hi, func):
	res = 0
	for i in range(lo, hi+1):
		for j in range(lo, hi+1):
			if j == i:
				continue
			else:
				for k in range(lo, hi+1):
					if k == i or k == j:
						continue
					else:
						res += func(i,j,k)
	return res

def QuadCLSum(lo, hi, func):
	res = 0
	for i in range(lo, hi+1):
		for j in range(lo, hi+1):
			if j == i:
				continue
			else:
				for k in range(lo, hi+1):
					if k == i or k == j:
						continue
					else:
						for l in range(lo, hi+1):
							if l == k or l== i or l == j:
								continue
							else:
								res += func(i,j,k,l)
	return res

def params(exp, cond, fit='both'):
	res = {'m1': 0,
			'm2': 0,
			'm3': 0,
			't12': 0,
			't13': 0,
			't23': 0,
			'd13': 0,
			'l': 0,
			'E': 0,
			'A': 0,
			'Ap': 0}
	if fit == "both":
		if exp == "nova":
			res['t12'] = np.arcsin(np.sqrt(0.31))
			res['t13'] = np.arcsin(np.sqrt(0.022))
			res['l'] = 810
			res['E'] = 1.9
			res['A'] = 7.56 * 10 ** (-5) * 2.75 * 1.9
			res['Ap'] = 3.815 * 10 ** (-5) * 2.75 * 1.9
			if cond == "IH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.36 * 10**(-3))
				res['m1'] = m31
				res['m2'] = m21 + m31
				res['m3'] = 0
				res['t23'] = np.arcsin(np.sqrt(0.56))
				res['d13'] = 285 / 360 * 2 * np.pi
			else:
				m21 = (7.39 * 10**(-5))
				m31 = (2.52 * 10**(-3))
				res['m1'] = 0
				res['m2'] = m21
				res['m3'] = m31
				res['t23'] = np.arcsin(np.sqrt(0.46))
				res['d13'] = 222 / 360 * 2 * np.pi
		if exp == "t2k":
			res['t12'] = np.arcsin(np.sqrt(0.31))
			res['t13'] = np.arcsin(np.sqrt(0.022))
			res['l'] = 295
			res['E'] = 0.6
			res['A'] = 7.56 * 10 ** (-5) * 3.4 * 0.6
			res['Ap'] = 3.815 * 10 ** (-5) * 3.4 * 0.6
			if cond == "IH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.46 * 10**(-3))
				res['m1'] = m31
				res['m2'] = m21 + m31
				res['m3'] = 0
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 285 / 360 * 2 * np.pi
			if cond == "NH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.56 * 10**(-3))
				res['m1'] = 0
				res['m2'] = m21
				res['m3'] = m31
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 222 / 360 * 2 * np.pi
	if fit == "t2k":
		if exp == "nova":
			res['t12'] = np.arcsin(np.sqrt(0.31))
			res['t13'] = np.arcsin(np.sqrt(0.022))
			res['l'] = 810
			res['E'] = 1.9
			res['A'] = 7.56 * 10 ** (-5) * 2.75 * 1.9
			res['Ap'] = 3.815 * 10 ** (-5) * 2.75 * 1.9
			if cond == "IH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.46 * 10**(-3))
				res['m1'] = m31
				res['m2'] = m21 + m31
				res['m3'] = 0
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 285 / 360 * 2 * np.pi
			if cond == "NH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.56 * 10**(-3))
				res['m1'] = 0
				res['m2'] = m21
				res['m3'] = m31
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 222 / 360 * 2 * np.pi
		if exp == "t2k":
			res['t12'] = np.arcsin(np.sqrt(0.31))
			res['t13'] = np.arcsin(np.sqrt(0.022))
			res['l'] = 295
			res['E'] = 0.6
			res['A'] = 7.56 * 10 ** (-5) * 3.4 * 0.6
			res['Ap'] = 3.815 * 10 ** (-5) * 3.4 * 0.6
			if cond == "IH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.46 * 10**(-3))
				res['m1'] = m31
				res['m2'] = m21 + m31
				res['m3'] = 0
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 285 / 360 * 2 * np.pi
			if cond == "NH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.56 * 10**(-3))
				res['m1'] = 0
				res['m2'] = m21
				res['m3'] = m31
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 222 / 360 * 2 * np.pi
	if fit == 'nova':
		if exp == "nova":
			res['t12'] = np.arcsin(np.sqrt(0.31))
			res['t13'] = np.arcsin(np.sqrt(0.022))
			res['l'] = 810
			res['E'] = 1.9
			res['A'] = 7.56 * 10 ** (-5) * 2.75 * 1.9
			res['Ap'] = 3.815 * 10 ** (-5) * 2.75 * 1.9
			if cond == "IH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.36 * 10**(-3))
				res['m1'] = m31
				res['m2'] = m21 + m31
				res['m3'] = 0
				res['t23'] = np.arcsin(np.sqrt(0.56))
				res['d13'] = 285 / 360 * 2 * np.pi
			else:
				m21 = (7.39 * 10**(-5))
				m31 = (2.52 * 10**(-3))
				res['m1'] = 0
				res['m2'] = m21
				res['m3'] = m31
				res['t23'] = np.arcsin(np.sqrt(0.46))
				res['d13'] = 222 / 360 * 2 * np.pi
		if exp == "t2k":
			res['t12'] = np.arcsin(np.sqrt(0.31))
			res['t13'] = np.arcsin(np.sqrt(0.022))
			res['l'] = 295
			res['E'] = 0.6
			res['A'] = 7.56 * 10 ** (-5) * 3.4 * 0.6
			res['Ap'] = 3.815 * 10 ** (-5) * 3.4 * 0.6
			if cond == "IH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.46 * 10**(-3))
				res['m1'] = m31
				res['m2'] = m21 + m31
				res['m3'] = 0
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 285 / 360 * 2 * np.pi
			if cond == "NH":
				m21 = (7.39 * 10**(-5))
				m31 = (2.56 * 10**(-3))
				res['m1'] = 0
				res['m2'] = m21
				res['m3'] = m31
				res['t23'] = np.arcsin(np.sqrt(0.55))
				res['d13'] = 222 / 360 * 2 * np.pi
	return res
			

def arbiangles(t12, t13, t23):
	return([t12, t13, t23])

def arbisangles(t14, t24, t34):
	return([t14, t24, t34])

def arbicp(d13, d14, d34):
	return([d13, d14, d34])

def arbimass(m1, m2, m3):
	return([m1, m2, m3])

def arbidensity(a, ap):
	return([a, ap])

def novaIHangles():
	return arbiangles(np.arcsin(np.sqrt(0.31))
		,np.arcsin(np.sqrt(0.022))
		,np.arcsin(np.sqrt(0.56)))

def novaNHangles():
	return arbiangles(np.arcsin(np.sqrt(0.31))
		,np.arcsin(np.sqrt(0.022))
		,np.arcsin(np.sqrt(0.46)))

def novaNHmass():
	return arbimass(0
		,np.sqrt(7.39 * 10 ** (-5))
		,np.sqrt(2.52 * 10 ** (-3)))

def novaIHmass():
	return arbimass(np.sqrt(2.36 * 10 ** (-3))
		,np.sqrt(7.39 * 10 ** (-5)) + np.sqrt(2.36 * 10 ** (-3))
		,0)

def novadensity():
	return([7.56 * 10 ** (-5) * 2.75 * 1.9
		,3.815 * 10 ** (-5) * 2.75 * 1.9])

def novaIHcp(d14, d34):
	return([285/360*np.pi*2, d14, d34])

def novaNHcp(d14, d34):
	return([222/360*np.pi*2, d14, d34])

def t2kIHangles():
	return arbiangles(np.arcsin(np.sqrt(0.31))
		,np.arcsin(np.sqrt(0.022))
		,np.arcsin(np.sqrt(0.55)))

def novaparam():
	return([810, 1.9])

def t2kNHangles():
	return arbiangles(np.arcsin(np.sqrt(0.31))
		,np.arcsin(np.sqrt(0.022))
		,np.arcsin(np.sqrt(0.55)))

def t2kNHmass():
	return arbimass(0
		,np.sqrt(7.39 * 10 ** (-5))
		,np.sqrt(2.56 * 10 ** (-3)))

def t2kIHmass():
	return arbimass(np.sqrt(2.46 * 10 ** (-3))
		,np.sqrt(7.39 * 10 ** (-5)) + np.sqrt(2.46 * 10 ** (-3))
		,0)

def t2kdensity():
	return([7.56 * 10 ** (-5) * 3.4 * 0.6
		,3.815 * 10 ** (-5) * 3.4 * 0.6])

def t2kIHcp(d14, d34):
	return([285/360*np.pi*2, d14, d34])

def t2kNHcp(d14, d34):
	return([222/360*np.pi*2, d14, d34])

def t2kparam():
	return([295, 0.6])