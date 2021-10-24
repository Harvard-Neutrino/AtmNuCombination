import numpy as np
import nuSQUIDSpy as nsq
import nuflux

def AtmOsc(ipnu, CosZenith, Enu, theta13, theta23, dm32, MO, theta12=0.563942, dm21=7.65e-05):
	E_min = 1.0e-1
	E_max = 4.0e2
	E_nodes = 100
	energy_nodes = nsq.geomspace(E_min,E_max,E_nodes)

	cth_min = -1.0
	cth_max = 1.0
	cth_nodes = 40
	cz_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

	neutrino_flavors = 3
	interactions = 'False'

	nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

	nsq_atm.Set_MixingAngle(0, 1, theta12)
	nsq_atm.Set_MixingAngle(0, 2, theta13)
	nsq_atm.Set_MixingAngle(1, 2, theta23)
	nsq_atm.Set_SquareMassDifference(1, dm21)
    if MO==1:
        nuSQ.Set_SquareMassDifference(2,dm32)
    elif MO==-1:
        nuSQ.Set_SquareMassDifference(2,mh*dm32+dm21)

    AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
    flux = nuflux.makeFlux('IPhonda2014_sk_solmin')
    for ic,cth in enumerate(cz_nodes):
        for ie,E in enumerate(energy_nodes):
            nu_energy = E*units.GeV
            nu_cos_zenith = cth
            AtmInitialFlux[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
            AtmInitialFlux[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
            AtmInitialFlux[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
            AtmInitialFlux[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
            AtmInitialFlux[ic][ie][0][2] = 0 # flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau
            AtmInitialFlux[ic][ie][1][2] = 0 # flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar

	nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
	nsq_atm.EvolveState()

	osc = np.array([])

	for nu, cz, E in zip(ipnu, CosZenith, Enu):
		neuflavor = int(abs(nu) / 2) % 6
		if nu<0: neutype = 1
		else neutype = 0
		osc = np.append(osc, nsq_atm.EvalFlavor(neuflavor, cz, E * units.GeV, neutype))

	return osc