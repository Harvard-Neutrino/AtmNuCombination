import numpy as np
import pandas as pd
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux
from params import *



units = nsq.Const()

Time = 3.0*365*24*60*60
meter_to_cm_sq = 1e4
Radian = 180. / np.pi

theta12 = np.arcsin(np.sqrt(0.304))
theta13 = np.arcsin(np.sqrt(0.02221))
theta23 = np.arcsin(np.sqrt(0.570))
m21 = 7.42e-5
m31 = 2.517e-3


########################    DATA   #####################################

# Define path to file (you may need to change this to match your system)
input_file = "neutrino_mc.csv"

# Load the file using pandas
input_data = pd.read_csv(input_file)



#######################  USE DATA    #####################################


# Define some energy bins (used throughout this notebook)
energy_bins_fine = np.logspace(1., 2., num=21)
energy_bins_course = np.logspace(1., 2., num=11)

# Define masks to identify different neutrino flavors
nue_mask = (np.abs(input_data["pdg"]) == 12)
numu_mask = (np.abs(input_data["pdg"]) == 14)
nutau_mask = (np.abs(input_data["pdg"]) == 16)


# Define masks to identify different flavor/interaction combinations.
nc_mask = input_data["current_type"] == 0
cc_mask = input_data["current_type"] == 1
nue_cc_mask = nue_mask & cc_mask
numu_cc_mask = numu_mask & cc_mask
nutau_cc_mask = nutau_mask & cc_mask



##################    Flux     ###########################

flux = nuflux.makeFlux('IPhonda2014_spl_solmin')

interactions = False

E_min = 1.0*units.GeV
E_max = 1.0e3*units.GeV
E_nodes = 200
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3

nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

#Initialize the flux
AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
flux = nuflux.makeFlux('IPhonda2014_spl_solmin')
for ic,cth in enumerate(nsq_atm.GetCosthRange()):
    for ie,E in enumerate(nsq_atm.GetERange()):
        nu_energy = E/units.GeV
        nu_cos_zenith = cth
        AtmInitialFlux[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
        AtmInitialFlux[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
        AtmInitialFlux[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
        AtmInitialFlux[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
        AtmInitialFlux[ic][ie][0][2] = 0.  # nutau
        AtmInitialFlux[ic][ie][1][2] = 0.  # nutau bar

        #AtmInitialFlux[ic][ie][0][2] = flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau
        #AtmInitialFlux[ic][ie][1][2] = flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar



        

######################### Evolution flux throught the Earth    ######################################
# nsq_atm.Set_MixingAngle(0,1, 33.44 * Radian)
# nsq_atm.Set_MixingAngle(0,2, 8.57 * Radian)
# nsq_atm.Set_MixingAngle(1,2, 49.2 * Radian)

# nsq_atm.Set_SquareMassDifference(1,7.42e-5)
# nsq_atm.Set_SquareMassDifference(2,2.517e-3)
nsq_atm.Set_MixingAngle(0, 1, theta12)
nsq_atm.Set_MixingAngle(0, 2, theta13)
nsq_atm.Set_MixingAngle(1, 2, theta23)
nsq_atm.Set_SquareMassDifference(1, m21)
nsq_atm.Set_SquareMassDifference(2, m31)

nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
#nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
nsq_atm.EvolveState()



##################################### Event distribution    #####################################

Erec_min = 1
Erec_max = 1000
NErec = 30
erec = np.logspace(np.log10(Erec_min), np.log10(Erec_max), NErec+1, endpoint = True)

dlter = (np.log10(Erec_max) - np.log10(Erec_min))/NErec

Crec_min = -1
Crec_max = 1
Ncrec = 10
crec = np.linspace(Crec_min, Crec_max, Ncrec+1, endpoint = True)

dltcr = (Crec_max - Crec_min)/Ncrec


rate_weight = np.zeros_like(input_data["weight"])
EventsBF = np.zeros((2, NErec, Ncrec))
EvBFEne = np.zeros((2, NErec))
EvBFCzen = np.zeros((2, Ncrec))

for i in range(len(rate_weight)):

    if input_data["pdg"][i] > 0 :
        neutype = 0
    else:
        neutype = 1
                
    if input_data["pid"][i] == 0 :
        neuint = 0
    elif input_data["pid"][i] == 1 :
        neuint = 1

    if np.abs(input_data["pdg"][i]) == 12:
        neuflavor = 0
    elif np.abs(input_data["pdg"][i]) == 14:
        neuflavor = 1
    elif np.abs(input_data["pdg"][i]) == 16:
        neuflavor = 2


    if input_data["reco_energy"][i] * units.GeV > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:
        rate_weight[i] = np.histogram(input_data["reco_energy"][cascade_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][cascade_mask])
        # binr = int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter)
        # binc = int((np.cos(input_data["reco_zenith"][i]) - Crec_min)/dltcr)

        # EventsBF[neuint][binr][binc] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

        
        # EvBFEne[neuint][binr] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

        # EvBFCzen[neuint][binc] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
input_data["rate_weight"] = rate_weight
energy_hist_truth, energy_bins_truth = np.histogram(input_data["reco_energy"][cascade_mask], bins = energy_bins_fine, weights = input_data["rate_weight"][cascade_mask])
print("get_rated_weight_truth: energy bins: \n", energy_bins_truth)
print("get_rated_weight_truth: energy rates: \n", energy_hist_truth)
# for e in range(NErec) :
#     print(erec[e], end="    ")
#     print(EvBFEne[0][e], end="   ")
#     print(EvBFEne[1][e], end="\n")
#     print(erec[e+1], end="    ")
#     print(EvBFEne[0][e], end="   ")
#     print(EvBFEne[1][e], end="\n")

# exit(0)


