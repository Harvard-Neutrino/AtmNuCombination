import numpy as np
import random
from math import pi, sqrt
import applications as ap
from particle import Particle


class nonFCReco:
	MeVtoGeV = 0.001
	def __init__(self, pdg, itype, distros, p, P):
		self.pdg = pdg
		self.distros = distros
		self.itype = itype
		self.TrueDirection = p/P
		self.TrueMomentum = P

		self.ReconDirection()
		self.ReconMomentum()

	def ReconDirection(self):
		ang=10
		# if self.itype == 14 or self.itype == 15: # PC-Stop, PC-Thru
		if self.itype == 14: # PC-Stop, PC-Thru
			ang = self.distros.Random('ang_pc')
		# elif self.itype == 16: # UpMu-Stop
		elif self.itype == 16 or self.itype == 15: # UpMu-Stop
			ang = self.distros.Random('ang_upmustop')
		elif self.itype == 17 or self.itype == 18: # UpMu-Thru
			ang = self.distros.Random('ang_upmuthru')

		ang=ang*pi/180.
		# Rodrigues way
		u = ap.RndVector()
		self.Direction = ap.RodRot(self.TrueDirection,u,ang)

	def ReconMomentum(self):
		reso = 0.01*(1.7+0.7/sqrt(self.TrueMomentum))
		bias = 0.0025 * self.TrueMomentum + 0.0015
		mHypothesis = Particle.from_pdgid(13).mass * self.MeVtoGeV
		p = (1-random.gauss(bias,reso)) * self.TrueMomentum
		pCorrSq = p**2+(Particle.from_pdgid(self.pdg).mass * self.MeVtoGeV)**2-mHypothesis**2
		if pCorrSq < 0:
			self.Momentum = 0
		else:
			self.Momentum = sqrt(pCorrSq)