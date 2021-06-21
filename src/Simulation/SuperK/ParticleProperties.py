import numpy as np
from math import sqrt, asin, atan, sin, cos
from applications import RodRot, RndVector, CrossProduct
from particle import Particle
import phasespace as ps

def MesonProduction(ip,mode):
    flag = 0
    if mode<30:
        modes = np.array([11,13,15,16,18,19,20,21,22,23,26])
    elif mode>30 and ip==2:
        modes = np.array([33,34,41,42,43,44,45,46])
    else:
        modes = np.array([41,44,45,46])
    for i in modes:
        if i==mode:
            flag=1

    return flag


# Explicitly store the primary lepton PDG ID
def LeptonPDG(NC,nu_flavour):
    l_id=0

    if NC:
        l_id=nu_flavour #No FCNC
    else:
        l_id = (abs(nu_flavour) - 1) * np.sign(nu_flavour)

    return l_id


# # Quicly decaying particles (< 5 ns)
# def DecayOf(PDG,E,p):

# 	return 0.



# # Pi0 decay
# def oldPi0Decay(e,p,dirv):
#     mpi0 = Particle.from_pdgid(111).mass * 0.001
#     sin_alpha_2 = 2
#     while abs(sin_alpha_2)>1:
#         E1 = np.random.rand() * e/2
#         E2 = e - E1
#         sin_alpha_2 = sqrt(mpi0**2 / (4*E1*E2)) # angle between photons
    
#     alpha = 2 * asin(sin_alpha_2)

#     alpha1 = atan(sin(alpha) / (cos(alpha) - E1/E2))
#     alpha2 = alpha - alpha1

#     u = np.array([0.,0.,0.])
#     u = RndVector()
#     dir1 = RodRot(dirv,u,alpha1)
#     dir2 = RodRot(dirv,u,-alpha2)
#     dir1 = dir1 / np.linalg.norm(dir1)
#     dir2 = dir2 / np.linalg.norm(dir2)

#     Eboth = np.array([E1,E2])
#     dirboth = np.array([dir1,dir2])

#     return Eboth, dirboth

def Pi0Decay(e,p,dirv):
    print('There we go!')
    MassPi0 = Particle.from_pdgid(111).mass
    MassGamma = 0.

    Momemtum4 = np.append(dirv * p, e) * 1000.

    decay = ps.nbody_decay(MassPi0, [MassGamma, MassGamma])
    weights, particles = decay.generate(n_events=1, boost_to=Momemtum4)

    p0 = particles['p_0'][:].numpy()
    p1 = particles['p_1'][:].numpy()

    dirboth = np.array([p0[0,0:3] / p0[0,3], p1[0,0:3] / p1[0,3]])
    Eboth = np.array([p0[0,3],p1[0,3]])

    return Eboth, dirboth


def InvariantMass(p1, p2, u1, u2):
    cth = np.dot(u1, u2)
    if cth>1: cth=1
    elif cth<-1: cth=-1
    return sqrt(2 * p1 * p2 * (1 - cth))

def ChargedKaonDecay():
    pass