import numpy as np
import pandas as pd
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux



units = nsq.Const()

Time = 3.0*365*24*60*60
meter_to_cm_sq = 1e4
Radian = 180. / np.pi


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
nsq_atm.Set_MixingAngle(0,1, 33.44 * Radian)
nsq_atm.Set_MixingAngle(0,2, 8.57 * Radian)
#nsq_atm.Set_MixingAngle(1,2, 0.8587)
nsq_atm.Set_MixingAngle(1,2, 49.2 * Radian)

nsq_atm.Set_SquareMassDifference(1,7.42e-5)
#nsq_atm.Set_SquareMassDifference(2,2.517e-3)
nsq_atm.Set_SquareMassDifference(2,2.517e-3)

nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
#nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
nsq_atm.EvolveState()



##################################### Event distribution    #####################################

Erec_min = 1
Erec_max = 1000000
NErec = 60
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

        binr = int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter)
        binc = int((np.cos(input_data["reco_zenith"][i]) - Crec_min)/dltcr)

        EventsBF[neuint][binr][binc] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

        
        EvBFEne[neuint][binr] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

        EvBFCzen[neuint][binc] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq



for e in range(NErec) :
    print(erec[e], end="    ")
    print(EvBFEne[0][e], end="   ")
    print(EvBFEne[1][e], end="\n")
    print(erec[e+1], end="    ")
    print(EvBFEne[0][e], end="   ")
    print(EvBFEne[1][e], end="\n")

exit(0)

#for c in range(Ncrec) :
#    print(crec[c], end="   ")
#    print(EvBFCzen[0][c], end="   ")
#    print(EvBFCzen[1][c], end="\n")
#    print(crec[c+1], end="   ")
#    print(EvBFCzen[0][c], end="   ")
#    print(EvBFCzen[1][c], end="\n")



#
#    print(crec[c+1], end="   ")
#    print(EvBFCzen[0][c], end="   ")
#    print(EvBFCzen[1][c], end="\n")
#



exit(0)



NT23 = 40
SqT23_min = 0.3
SqT23_max = 0.7
SqT23 = np.linspace(SqT23_min, SqT23_max, NT23+1, endpoint = True)
T23 = np.arcsin(np.sqrt(SqT23))


NDM31 = 40
DM31_max = 2.8e-3
DM31_min = 2.3e-3
DM31 = np.linspace(DM31_min, DM31_max, NDM31+1, endpoint = True)


Events = np.zeros((2, NDM31, NT23, NErec, Ncrec))
EventsMass = np.zeros((2, NDM31, NErec, Ncrec))
#Events = np.zeros((2, NT23, NErec, Ncrec))
#Events = np.zeros((2, NErec, Ncrec))

#for m in range(NDM31):
#    for t in range(NT23):




#for t in range(NT23):
for m in range(NDM31):

    #for t in range(NT23):

    ######################### Evolution flux throught the Earth    ######################################
    nsq_atm.Set_MixingAngle(0,1, 33.44 * Radian)
    nsq_atm.Set_MixingAngle(0,2, 0.1496)
    nsq_atm.Set_MixingAngle(1,2, 0.865743)
    nsq_atm.Set_SquareMassDifference(1,7.42e-5)
    nsq_atm.Set_SquareMassDifference(2,DM31[m])


    nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
    #nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
    nsq_atm.EvolveState()



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


        if input_data["reco_energy"][i] > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:
            
            binr = int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter)
            binc = int((np.cos(input_data["reco_zenith"][i]) - Crec_min)/dltcr)

            EventsMass[neuint][m][binr][binc] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

            #Events[neuint][m][t][binr][binc] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

        


Xi = np.zeros((NDM31))

for m in range(NDM31):

    a = 0
    for s in range(2):
        for i in range(NErec) : 
            for j in range(Ncrec) : 
                
                if EventsBF[s][i][j] > 0 :
                    Xi[m] += (EventsMass[s][m][i][j] - EventsBF[s][i][j]) * (EventsMass[s][m][i][j] - EventsBF[s][i][j]) / EventsBF[s][i][j]
                    a += (EventsMass[s][m][i][j] - EventsBF[s][i][j]) * (EventsMass[s][m][i][j] - EventsBF[s][i][j]) / EventsBF[s][i][j]


    print(DM31[m], a)

for m in range(NDM31):
    print(DM31[m], Xi[m])




#Xi = np.zeros((NT23, NDM31))
#
#for m in range(NDM31):
#    for t in range(NT23):
#        a = 0
#        for s in range(2):
#            for i in range(NErec) : 
#                for j in range(Ncrec) : 
#                
#                    if EventsBF[s][i][j] > 0 :
#                        Xi[t][m] += (Events[s][m][t][i][j] - EventsBF[s][i][j]) * (Events[s][m][t][i][j] - EventsBF[s][i][j]) / EventsBF[s][i][j]
#                        a += (Events[s][m][t][i][j] - EventsBF[s][i][j]) * (Events[s][m][t][i][j] - EventsBF[s][i][j]) / EventsBF[s][i][j]
#
#
#        print(DM31[m], SqT23[t], a)
#
#for m in range(NDM31):
#    for t in range(NT23):
#        print(SqT23[t], DM31[m], Xi[t][m])


#for e in range(NErec) :
#    for c in range(Ncrec) :
#        print(erec[e], end="    ")
#        print(crec[c], end="   ")
#        print(EventsBF[0][e][c], end="   ")
#        print(EventsBF[1][e][c], end="   ")
#        print(Events[0][0][0][e][c], end="   ")
#        print(Events[1][0][0][e][c], end="\n")
#
#        print(erec[e], end="    ")
#        print(crec[c+1], end="   ")
#        print(EventsBF[0][e][c], end="   ")
#        print(EventsBF[1][e][c], end="   ")
#        print(Events[0][0][0][e][c], end="   ")
#        print(Events[1][0][0][e][c], end="\n")
#
#exit(0)

#Xi = np.zeros((NDM31, NT23))

#for m in range(NDM31):
#    for t in range(NT23):
#        
#        for s in range(2):
#            for i in range(NErec) : 
#                for j in range(Ncrec) : 
#                    if EventsBF[s][i][j] > 0 :
#                        Xi[m][t] += (Events[s][m][t][i][j] - EventsBF[s][i][j]) * (Events[s][m][t][i][j] - EventsBF[s][i][j]) / EventsBF[s][i][j]



#for m in range(NDM31):
#    for t in range(NT23):
#        print(DM31[m], SqT23[t], Xi[m][t])
