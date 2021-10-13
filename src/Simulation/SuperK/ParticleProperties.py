import numpy as np
from math import sqrt, asin, atan, sin, cos
from applications import RodRot, RndVector, CrossProduct
from particle import Particle

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



def InvariantMass(p1, p2, u1, u2):
    cth = np.dot(u1, u2)
    if cth>1: cth=1
    elif cth<-1: cth=-1
    return sqrt(2 * p1 * p2 * (1 - cth))
