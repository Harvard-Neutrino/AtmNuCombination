import numpy as np
import ParticleProperties as pp
from math import sqrt, asin, atan, sin, cos, pi


def ClearBadHadrons(nhad, norm, E, pdg, dire, mom, p):
	boo = np.where(norm!=norm)
	bad_hads = np.size(boo)

	new_norm = np.delete(norm,boo)
	new_E = np.delete(E,boo)
	new_p = np.delete(p,boo)
	new_pdg = np.delete(pdg,boo)

	new_dire = np.delete(dire,boo,0)
	new_mom = np.delete(mom,boo,0)

	new_nhad = nhad - bad_hads

	return new_nhad, new_norm, new_E, new_pdg, new_dire, new_mom, new_p


def RodRot(v,u,ang):
    vrot = np.array([0.,0.,0.])
    vrot[0] = (cos(ang)+u[0]**2*(1-cos(ang)))*v[0] + (u[0]*u[1]*(1-cos(ang))-u[2]*sin(ang))*v[1] + (u[0]*u[2]*(1-cos(ang)+u[1]*sin(ang)))*v[2]
    vrot[1] = (cos(ang)+u[1]**2*(1-cos(ang)))*v[1] + (u[0]*u[1]*(1-cos(ang))+u[2]*sin(ang))*v[0] + (u[1]*u[2]*(1-cos(ang)-u[0]*sin(ang)))*v[2]
    vrot[2] = (cos(ang)+u[2]**2*(1-cos(ang)))*v[2] + (u[0]*u[2]*(1-cos(ang))-u[1]*sin(ang))*v[0] + (u[1]*u[2]*(1-cos(ang)+u[0]*sin(ang)))*v[1]
    vrot = vrot / np.linalg.norm(vrot)
    return vrot

# Random unit vector
def RndVector():
    ph = 2*pi*np.random.rand()
    th = pi*np.random.rand()
    u = np.array([cos(ph)*sin(th),sin(ph)*sin(th),cos(th)])
    return u

def CrossProduct(a,b):
    c = np.array([0.,0.,0.])
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    modc = np.linalg.norm(c)
    c = c / modc
    return modc, c

def Azimuth(a, b):
    return asin(a/(sqrt(1-b**2)))