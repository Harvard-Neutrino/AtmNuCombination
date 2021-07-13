import numpy as np
from particle import Particle
import ParticleProperties as pp


def SKIP(pdg):
	if abs(pdg)==11:
		return 2
	elif pdg==22:
		return 2
	elif abs(pdg)==13:
		return 3
	elif abs(pdg)==211:
		return 3
	else:
		return 3

def TrueRingConstructor(pdg, E, P, Pv):

	nring=0
	ring_pdg = np.array([])
	ring_ip = np.array([])
	ring_energy = np.array([])
	ring_momentum = np.array([])
	ring_direction = np.array([])

	# Loop over event final state particles checking if they're capable of producing a ring
	for i, p in enumerate(P):
		if p>0.001:
			if (pdg[i]==22 or (Particle.from_pdgid(pdg[i]).charge!=0) and pdg[i]!=-211):

				Ethr = (E[i] + 0.03) / 1.33
			
				if p > Ethr:
					ring_pdg = np.append(ring_pdg, pdg[i])
					ring_ip = np.append(ring_ip, SKIP(pdg[i]))
					# print(pdg[i], '--', SKIP(pdg[i]))
					ring_energy = np.append(ring_energy, E[i])
					ring_momentum = np.append(ring_momentum, p)

					if not nring:
						ring_direction = Pv[i] / p
					else:
						ring_direction = np.vstack((ring_direction, Pv[i] / p))

					nring += 1

	return nring, ring_pdg, ring_ip, ring_energy, ring_momentum, ring_direction

def numberOfNeutrons(pdgs):
	nn = 0
	for i in pdgs:
		if i==2112:
			nn += 1
	return nn