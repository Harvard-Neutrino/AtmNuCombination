import numpy as np
# from particle import Particle
# import ParticleProperties as pp
# import math
from math import log10, pi, sqrt, asin
import ROOT as root
# import mulself.TrueIProcessing
# import random
# import applications as ap
# from Apps import TrueRingConstructor

class Neutrino:

	deg2rad = pi/180.

	def __init__(self, htree, ipnu, E, P_v, mode, cc, coh, dis, nc):
		self.Mode = mode
		self.Energy = E
		self.Flavour = ipnu
		self.Direction = P_v / E
		self.CosZenith = self.Direction[2]
		self.Azi = asin(self.Direction[1]/(sqrt(1-self.Direction[2]**2)))
		self.Momenta = P_v
		self.CC = cc
		self.COH = coh
		self.DIS = dis
		self.NC = nc
		self.htree = htree

	def LeptonPDG(self):
		if self.NC :
			return self.Flavour #No FCNC
		else:
			return (abs(self.Flavour) - 1) * np.sign(self.Flavour)

	def Flux(self): # custom ROOT flux tree based on HKKM's flux tables
		CZ_points = 10
		E_points = 20
		Azi_step = 30*self.deg2rad

		E_per_Azi = 101
		Azi_per_CZ = 12

		E_index = int(E_points*(log10(self.Energy)+1))
		CZ_index = int(CZ_points*(1-self.CosZenith)) # decreasing order in cos zenith
		Azi_index = int(self.Azi / Azi_step)

		tree_index0 = E_index + E_per_Azi*Azi_index + E_per_Azi*Azi_per_CZ*CZ_index
		tree_index1 = tree_index0 + 1

		tree = self.htree
		tree.GetEntry(tree_index0)
		E0 = tree.logE
		flux0 = np.array([tree.flux[0],tree.flux[1],tree.flux[2],tree.flux[3]])
		tree.GetEntry(tree_index1)
		E1 = tree.logE
		flux1 = np.array([tree.flux[0],tree.flux[1],tree.flux[2],tree.flux[3]])

		flux = np.zeros(4)
		flux[0] = np.interp(log10(self.Energy), [E0,E1], [flux0[0],flux1[0]])
		flux[1] = np.interp(log10(self.Energy), [E0,E1], [flux0[1],flux1[1]])
		flux[2] = np.interp(log10(self.Energy), [E0,E1], [flux0[2],flux1[2]])
		flux[3] = np.interp(log10(self.Energy), [E0,E1], [flux0[3],flux1[3]])

		return flux	

