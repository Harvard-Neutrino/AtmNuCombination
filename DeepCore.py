import numpy as np
import pandas as pd
#import matplotlib
#import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq
import nuflux



units = nsq.Const()

Time = 5.0*365*24*60*60
meter_to_cm_sq = 1e4


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

flux = nuflux.makeFlux('honda2006')

interactions = False

E_min = 10.0*units.GeV
E_max = 1.0e3*units.PeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3

nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

#Initialize the flux
AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
flux = nuflux.makeFlux('honda2006')
for ic,cth in enumerate(nsq_atm.GetCosthRange()):
    for ie,E in enumerate(nsq_atm.GetERange()):
        nu_energy = E/units.GeV
        nu_cos_zenith = cth
        AtmInitialFlux[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
        AtmInitialFlux[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
        AtmInitialFlux[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
        AtmInitialFlux[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
        AtmInitialFlux[ic][ie][0][2] = flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau
        AtmInitialFlux[ic][ie][1][2] = flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar



#print('Before evolution')
#for ic,cth in enumerate(nsq_atm.GetCosthRange()):
#for ie,E in enumerate(nsq_atm.GetERange()):
#    print("{:.2e}".format(E), "{:.2e}".format(AtmInitialFlux[0][ie][0][1]))
    #print(cth[0], "{:.2e}".format(E), "{:.2e}".format(AtmInitialFlux[0][ie][0][1]))
        


######################### Evolution flux throught the Earth    ######################################
nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
#nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
nsq_atm.EvolveState()

#print( "{:.2e}".format(nsq_atm.EvalFlavor(0, -1., 1e10, 0)))

#print('After evolution')
#for ic,cth in enumerate(nsq_atm.GetCosthRange()):
#for ie,E in enumerate(nsq_atm.GetERange()):
#    print("{:.2e}".format(E), "{:.2e}".format(nsq_atm.EvalFlavor(1, -1., E, 0)))



##################################### Event distribution    #####################################

Erec_min = 1
Erec_max = 1000000
NErec = 60
erec = np.logspace(np.log10(Erec_min), np.log10(Erec_max), NErec+1, endpoint = True)


dlter = (np.log10(Erec_max) - np.log10(Erec_min))/NErec



rate_weight = np.zeros_like(input_data["weight"])
inter = np.zeros((2, neutrino_flavors, NErec))

for i in range(len(rate_weight)):

    if input_data["pdg"][i] > 0 :
        neutype = 0
    else:
        neutype = 1
        
    if np.abs(input_data["pdg"][i]) == 12:
        neuflavor = 0
    elif np.abs(input_data["pdg"][i]) == 14:
        neuflavor = 1
    elif np.abs(input_data["pdg"][i]) == 16:
        neuflavor = 2


    if input_data["reco_energy"][i] > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:
        #binr = int( (np.log10(input_data["true_energy"][i]) - np.log10(Erec_min))/dlter)
        binr = int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter)

        #if input_data["pid"][i] == 0 :
        #    inter[0][neuflavor][binr] += 1
        #else :
        #    inter[1][neuflavor][binr] += 1


        if input_data["pid"][i] == 0 :
            inter[0][neuflavor][binr] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
        else :
            inter[1][neuflavor][binr] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq


    

    if input_data["true_energy"][i] * units.GeV < E_min or input_data["true_energy"][i] * units.GeV > E_max :
        rate_weight[i] = 0
        continue

    rate_weight[i] = input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq




reco_energy_bins = np.linspace(5,1000,100)

hist = np.histogram(input_data["reco_energy"], weights = rate_weight, bins = reco_energy_bins)

#print(hist)

for i in range(NErec) : 
    print(erec[i],  inter[0][0][i],  inter[0][1][i],  inter[0][2][i],  inter[1][0][i],  inter[1][1][i],  inter[1][2][i])
    print(erec[i+1],  inter[0][0][i],  inter[0][1][i],  inter[0][2][i],  inter[1][0][i],  inter[1][1][i],  inter[1][2][i])

#print("nu_e Cascade: ", inter[0][0])
#print("nu_m Cascade: ", inter[0][1])
#print("nu_t Cascade: ", inter[0][2])
#print("nu_e Track: ", inter[1][0])
#print("nu_m Track: ", inter[1][1])
#print("nu_t Track: ", inter[1][2])

