import h5py
import numpy as np
import math
import AtmOscillator as atm
import nuSQuIDS as nsq
import nuSQUIDSTools
import argparse


def FluxFactor(i, flavor, nue, nueb, numu, numub):
    j = int(abs(flavor) / 2) % 6
    factor = 1.0
    if i==j or i==1 and j==2:
        factor = 1.0
    elif i==0 and j>=1:
        if flavor>0:
            factor = nue / numu
        else:
            factor = nueb / numub
    elif i==1 and j==0:
        if flavor>0:
            factor = numu / nue
        else:
            factor = numub / nueb
    
    return factor



parser = argparse.ArgumentParser()
parser.add_argument("unosc_hdf5filename", type=str, nargs='?', default='../Simulation/SuperK/data/testfcmc.hdf5')
args = parser.parse_args()
filename = args.unosc_hdf5filename

with h5py.File(filename, 'r') as hf:
    print('Opening unoscillated file...')
    Enu = np.array(hf['pnu'][()]) # true neutrino energy
    neu = np.array(hf['ipnu'][()]) # true neutrino flavour
    dirnuX = np.array(hf['dirnuX'][()]) # true neutrino X-direction
    dirnuY = np.array(hf['dirnuY'][()]) # true neutrino Y-direction
    dirnuZ = np.array(hf['dirnuZ'][()]) # true neutrino Z-direction
    
    mode = np.array(hf['mode'][()]) # neutrino interaction mode
    flx_nue = np.array(hf['fluxho_nue'][()]) # nue flux
    flx_nueb = np.array(hf['fluxho_nueb'][()]) # nuebar flux
    flx_numu = np.array(hf['fluxho_numu'][()]) # numu flux
    flx_numub = np.array(hf['fluxho_numub'][()]) # numubar flux

    evis = np.array(hf['evis'][()]) # reconstructed neutrino energy
    recodirZ = np.array(hf['recodirZ'][()]) # reconstructed neutrino Z-direction
    itype = np.array(hf['itype'][()]) # reconstructed event type


''' Oscillation parameters grid, to xml/card file? '''
# Oscillation parameters
#Usually fixed
th12 = 0.563942
dm21 = 7.65e-05

#Usually not fixed
th13_0 = math.asin(math.sqrt(0.018))
th13_end = math.asin(math.sqrt(0.018))
nth13 = 1
th13_step = (th13_end-th13_0)/(float(nth13))

th23_0 = math.asin(math.sqrt(0.35))
th23_end = math.asin(math.sqrt(0.75))
nth23 = 20
th23_step = (th23_end-th23_0)/(float(nth23))

dm32_0 = 0.0015
dm32_end = 0.004
ndm32 = 25
dm32_step = (dm32_end-dm32_0)/(float(ndm32))

dcp_0  = 3*math.pi/2
dcp_end = 3*math.pi/2
ndcp = 1
dcp_step = (dcp_end-dcp_0)/(float(ndcp))

mh_0   = 1
mh_end = -1
nmh = 2
mh_step = (mh_end-mh_0) / (float(nmh)-1)



osc_points = nmh*ndcp*ndm32*nth13*nth23

# outfile to write
# with h5py.File('osc_hdf5', 'w') as hf:

w = np.zeros(osc_points).reshape(nmh,ndcp,ndm32,nth13,nth23)

units = nsq.Const()

for k,(nu,E,cz,mod) in enumerate(zip(neu, Enu, dirnuZ, mode)):
    # Get P_{x->ipnu} probabilities
    weight = 0.0
    if nu>0:
        nuSQ = nsq.nuSQUIDS(3,nsq.NeutrinoType.neutrino)
    elif nu<0:
        nuSQ = nsq.nuSQUIDS(3,nsq.NeutrinoType.antineutrino)
    else:
        print('What?! No identified neutrino flavour')
    nuSQ.Set_E(E*units.GeV)
    #print('*** Earth Atmospheric Neutrino Osc ***')
    nuSQ.Set_Body(nsq.EarthAtm())
    zenith = np.arccos(cz)
    nuSQ.Set_Track(nsq.EarthAtm().Track(zenith))
    nuSQ.Set_rel_error(1.0e-4);
    nuSQ.Set_abs_error(1.0e-4);

    for j in range(nmh):
        mh=mh_0+j*mh_step
        for k in range(ndcp):
            dcp=dcp_0+k*dcp_step
            for l in range(ndm32):
                dm32=dm32_0+l*dm32_step
                for m in range(nth13):
                    th13=th13_0+m*th13_step
                    for n in range(nth23):
                        th23=th23_0+n*th23_step
                        weight = 0.0
                        if abs(mod) < 30:
                            for i in range(2):
                                in_state = np.zeros(3)
                                in_state[i] = 1
                                Ffactor = FluxFactor(i, nu, flx_nue[k], flx_nueb[k], flx_numu[k], flx_numub[k])
                                prob = atm.SMOsc(nuSQ,nu,in_state,th13,th23,th12,dm21,dm32,dcp,mh)
                                weight += prob*Ffactor
                        else:
                            weight = 1.0
                        w[j][k][l][m][n] = weight

    # print(w)



