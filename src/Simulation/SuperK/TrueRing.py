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

def TrueRingConstructor(PDGlep, PDGhad, Elep, Ehad, Plep, Phad, Plepv, Phadv):
	pdg = np.append(PDGlep, PDGhad)
	E = np.append(Elep, Ehad)
	P = np.append(Plep, Phad)
	Pv = np.vstack((Plepv, Phadv))

	nring=0
	ring_pdg = np.array([0.])
	ring_ip = np.array([0.])
	ring_energy = np.array([0.])
	ring_momentum = np.array([0.])
	ring_direction = np.array([0.,0.,0.])

	# Loop over event final state particles checking if they're capable of producing a ring
	for i, p in enumerate(P):
		if p>0.001:
			if (pdg[i]==22 or (Particle.from_pdgid(pdg[i]).charge!=0 and abs(pdg[i])!=15)):

				Ethr = (E[i] + 0.03) / 1.33
			
				if p > Ethr:
					if not nring:
						ring_pdg[0] = pdg[i]
						ring_ip[0] = SKIP(pdg[i])
						ring_energy[0] = E[i]
						ring_momentum[0] = p
						ring_direction = Pv[i] / p
					else:
						ring_pdg = np.append(ring_pdg, pdg[i])
						ring_ip = np.append(ring_ip, SKIP(pdg[i]))
						ring_energy = np.append(ring_energy, E[i])
						ring_momentum = np.append(ring_momentum, p)
						ring_direction = np.vstack((ring_direction, Pv[i] / p))
					nring += 1

			# Example for introducing fast (O(ns)) particle decays
			elif pdg[i]==111: 

				E_pi0_gammas, dir_pi0_gammas = pp.Pi0Decay(E[i], p, Pv[i]/p)
				pdg_pi0_gammas = np.array([22,22])

				for j, gamma in enumerate(pdg_pi0_gammas):
					if E_pi0_gammas[j] > 0.03:
						if not nring:
							ring_pdg[0] = gamma
							ring_ip[0] = SKIP(gamma)
							ring_energy[0] = E_pi0_gammas[j]
							ring_momentum[0] = E_pi0_gammas[j]
							ring_direction = dir_pi0_gammas[j]
						else:
							ring_pdg = np.append(ring_pdg, gamma)
							ring_ip = np.append(ring_ip, SKIP(gamma))
							ring_energy = np.append(ring_energy, E_pi0_gammas[j])
							ring_momentum = np.append(ring_momentum, E_pi0_gammas[j])
							ring_direction = np.vstack((ring_direction,dir_pi0_gammas[j]))
						nring += 1

	return nring, ring_pdg, ring_ip, ring_energy, ring_momentum, ring_direction