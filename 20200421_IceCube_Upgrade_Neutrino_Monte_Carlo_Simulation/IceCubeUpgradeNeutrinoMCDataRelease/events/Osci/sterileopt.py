import math
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.linalg as la
import random

import sterileU
import sterileutil as util
import sterileprob

# define best fit values we want to optimize on
t2kNHfit = [0.0675, 0.03]
novaNHfit = [0.044, 0.0375]
t2kNHCI = [0.0085, 0.024]
novaNHCI = [0.0075, 0.0135]

# define probability functions with parameters we want opt on
def t2kNH(vt14, vt24, vd13, vd14):
	return [sterileprob.prob(2, 1, 't2k', 'IH', 't2k', m4 = 1, t14 = vt14, t24 = vt24,\
	 t34= 0, mand13 = True, vald13 = vd13,manE = False, valE = 0.9, d14 = vd14, d34 = 0, anti = False, avg = True), \
	sterileprob.prob(2, 1, 't2k', 'IH', 't2k', m4 = 1, t14 = vt14, t24 = vt24,\
	 t34= 0, mand13 = True, vald13 = vd13,manE = False, valE = 0.9, d14 = vd14, d34 = 0, anti = True, avg = True)]
def novaNH(vt14, vt24, vd13, vd14):
	return [sterileprob.prob(2, 1, 'nova', 'IH', 't2k', m4 = 1, t14 = vt14, t24 = vt24,\
	 t34= 0, mand13 = True, vald13 = vd13,manE = False, valE = 0.9, d14 = vd14, d34 = 0, anti = False, avg = True),\
	sterileprob.prob(2, 1, 'nova', 'IH', 't2k', m4 = 1, t14 = vt14, t24 = vt24,\
	 t34= 0, mand13 = True, vald13 = vd13,manE = False, valE = 0.9, d14 = vd14, d34 = 0, anti = True, avg = True)]

def NE(population = 7000, generation = 75, champion = 1, survivor = 59, \
		theta_mut = 0.01, delta_mut= 0.01):

	# define error function
	def error(vt14, vt24, vd13, vd14):
		t2kErr = ((t2kNH(vt14, vt24, vd13, vd14)[0] - t2kNHfit[0])**2 +\
						 (t2kNH(vt14, vt24, vd13, vd14)[1] - t2kNHfit[1])**2)
		novaErr = ((novaNH(vt14, vt24, vd13, vd14)[0] - novaNHfit[0])**2 +\
						 (novaNH(vt14, vt24, vd13, vd14)[1] - novaNHfit[1])**2)
		return t2kErr + novaErr 
	def errorCI(vt14, vt24, vd13, vd14):
		t2kErr = (np.absolute(t2kNH(vt14, vt24, vd13, vd14)[0] - t2kNHfit[0]) / t2kNHCI[0] \
				+ np.absolute(t2kNH(vt14, vt24, vd13, vd14)[1] - t2kNHfit[1]) / t2kNHCI[1])**2
		novaErr = (np.absolute(novaNH(vt14, vt24, vd13, vd14)[0] - novaNHfit[0]) / novaNHCI[0] \
				+ np.absolute(novaNH(vt14, vt24, vd13, vd14)[1] - novaNHfit[1]) / novaNHCI[1])**2
		return np.sqrt(t2kErr + novaErr)

	# define mutation function
	def mutate(parent):
		child = parent[:]
		min0 = child[0] - 0
		max0 = 0.2 - child[0]
		child[0] += np.random.uniform(-min(theta_mut, min0), min(theta_mut, max0))
		min1 = child[1] - 0
		max1 = 0.2 - child[1]
		child[1] += np.random.uniform(-min(theta_mut, min1), min(theta_mut, max1))
		child[2] += np.random.uniform(-delta_mut, delta_mut)
		child[3] += np.random.uniform(-delta_mut, delta_mut)
		return child

	# initialize entire populations
	pop = [[0, 0, 0, 0, 0]]
	for i in range(0, population - 1):
		pop.append([0, 0, 0, 0, 0])

	# randomize populations
	sample1 = np.random.uniform(0,0.1, (population,))
	sample2 = np.random.uniform(0,0.1, (population,))
	sample3 = np.random.uniform(-np.pi/2,np.pi/2, (population,))
	sample4 = np.random.uniform(-np.pi/2,np.pi/2, (population,))
	for i in range(0, len(pop)):
		pop[i][0] = sample1[i]
		pop[i][1] = sample2[i]
		pop[i][2] = sample3[i]
		pop[i][3] = sample4[i]

	

	# main loop of length generation
	for j in range(0, generation):
		print("\n in generation", j)
		for i in range(0, len(pop)):
			pop[i][4] = errorCI(pop[i][0], pop[i][1],\
							 pop[i][2], pop[i][3])
		pop.sort(key=lambda x:x[4])

		if j == generation - 1:
			print("FINISHED TRAINING WITH \n")
			print("champion is")
			print(pop[0])
			print("with fit values")
			print(t2kNH(pop[0][0], pop[0][1], pop[0][2], pop[0][3]))
			print(novaNH(pop[0][0], pop[0][1], pop[0][2], pop[0][3]))
			return

		# After sorting, select champ and survivors
		savesurv = []
		for champ in range(0, champion):
			savesurv.append(pop[champ])
		for surv in range(0, survivor):
			savesurv.append(pop[champion + surv])

		# first generate random picks from champion + survivors
		offspring = np.random.randint(champion + survivor,\
										size = population - champion)

		# first copy champions into next generation
		for champ in range(0, champion):
			pop[champ] = savesurv[champ][:]

		# now mutate the survivors according to the offspring index
		for child in range(0, population - champion):
			pop[child + champion] = mutate(savesurv[offspring[child]])

		print("\n the champion is now\n", pop[0])

		

# The following are used to document successful values

def error(vt14, vt24, vd13, vd14):
		t2kErr = (np.absolute(t2kNH(vt14, vt24, vd13, vd14)[0] - t2kNHfit[0])/ t2kNHCI[0] +\
						 np.absolute(t2kNH(vt14, vt24, vd13, vd14)[1] - t2kNHfit[1]) / t2kNHCI[1])
		novaErr = (np.absolute(novaNH(vt14, vt24, vd13, vd14)[0] - novaNHfit[0]) / novaNHCI[0] +\
						 np.absolute(novaNH(vt14, vt24, vd13, vd14)[1] - novaNHfit[1]) / novaNHCI[1])
		return t2kErr + novaErr

def errorCI(vt14, vt24, vd13, vd14):
		t2kErr = (np.absolute(t2kNH(vt14, vt24, vd13, vd14)[0] - t2kNHfit[0]) / t2kNHCI[0] \
				+ np.absolute(t2kNH(vt14, vt24, vd13, vd14)[1] - t2kNHfit[1]) / t2kNHCI[1])**2
		novaErr = (np.absolute(novaNH(vt14, vt24, vd13, vd14)[0] - novaNHfit[0]) / novaNHCI[0] \
				+ np.absolute(novaNH(vt14, vt24, vd13, vd14)[1] - novaNHfit[1]) / novaNHCI[1])**2
		return np.sqrt(t2kErr + novaErr)
def record(vt14, vt24, vd13, vd14):
	print("T2K", t2kNH(vt14, vt24, vd13, vd14))
	print("NOvA",novaNH(vt14, vt24, vd13, vd14))
	print("literature T2K:", t2kNHfit)
	print("literature NOvA:", novaNHfit)
	print("literature T2K 1CI:", t2kNHCI)
	print("literature NOvA 1CI:", novaNHCI)
	print("error,", error(vt14, vt24, vd13, vd14))

# These correspond to the older formula with opposite conjugates, different fit values
# record(0.22122758024394812, 0.11112785291907522, -0.24523177491810652, 1.4664307210456229)
# record(0.19183724533780158, 0.16709390128763782, -0.01547296343206932, 1.4640815753874934)
# record(0.1704502851872749, 0.18367733438823333, -0.06855384197330401, 1.5197067798275956)

# These are the new ones, deifferent fit values
# record(0.22122758024394812, 0.11112785291907522, 0.24523177491810652, -1.4664307210456229)
# record(0.19183724533780158, 0.16709390128763782, 0.01547296343206932, -1.4640815753874934)
# record(0.1704502851872749, 0.18367733438823333, 0.06855384197330401, -1.5197067798275956)
# record(0.09224084389302983, 0.31291614134920936, -0.03169622608432427, -1.4431026996218623) #error is 2.4214100826270746e-10

################################

# with t2k fit values these are different, these are the new solutions
# with t2k fit values but with sterile contribution averaged out
# record(0.2734963292488098, 0.3538208413398849, -1.70004508832868, -0.05162353132050409)
# without averaging out
# record(0.1782208531875613, 0.19173996859458908, -0.0663900618275011, 1.7562631105700481)


# now the same but with IH
# sterile averaged out
# not averaged out

# NE()
record(0.18702364582877923, 0.18166313457529676, -1.467247237927259, -0.5238424634332328)