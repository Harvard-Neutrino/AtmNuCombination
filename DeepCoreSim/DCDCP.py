import numpy as np
import pandas as pd
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux
import sys
from scipy.optimize import minimize

units = nsq.Const()

Time = 3.0*365*24*60*60
meter_to_cm_sq = 1e4
Radian = np.pi / 180.


########################    READ DATA   #####################################

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


############################### Oscillation Parameters ###############################
NT23 = 24
SqT23_min = 0.305
SqT23_max = 0.705
SqT23 = np.linspace(SqT23_min, SqT23_max, NT23+1, endpoint = True)
SqT23BF = SqT23[17]
T23 = np.arcsin(np.sqrt(SqT23))
T23BF = T23[16]


NDM31 = 24
DM31_max = 3.0e-3
DM31_min = 2.0e-3
DM31 = np.linspace(DM31_min, DM31_max, NDM31+1, endpoint = True)
DM31BF = DM31[12]

NDCP = 19
DCP_min = 0
DCP_max = 2*np.pi
DCP = np.linspace(DCP_min, DCP_max, NDCP+1, endpoint = True)
DCPBF = DCP[13]



########################### Selected bin ###############################################

Bin = int(sys.argv[1]) + 9000 
BT23 = int(Bin/(NDM31 * NDCP))
BM31 = int((Bin - BT23 * NDM31 * NDCP)/NDCP)
BDCP = int(Bin - BT23 * NDM31 * NDCP - BM31 * NDCP)


##################    Flux     ###########################

flux = nuflux.makeFlux('IPhonda2014_spl_solmin')

interactions = False

E_min = 1.0*units.GeV
E_max = 1.0e3*units.GeV
E_nodes = 100
energy_nodes = nsq.logspace(E_min,E_max,E_nodes)

cth_min = -1.0
cth_max = 1.0
cth_nodes = 40
cth_nodes = nsq.linspace(cth_min,cth_max,cth_nodes)

neutrino_flavors = 3

##Define the nuSquids function 
nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

#Initialize the flux
AtmInitialFlux = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
AtmIninue = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
AtmIninueBar = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
AtmIninum = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
AtmIninumBar = np.zeros((len(cth_nodes),len(energy_nodes),2,neutrino_flavors))
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


        AtmIninue[ic][ie][0][0] = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
        
        AtmIninueBar[ic][ie][1][0] = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue bar
        
        AtmIninum[ic][ie][0][1] = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
        
        AtmIninumBar[ic][ie][1][1] = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
        
        


        

######################### Evolution flux throught the Earth    ######################################
nsq_atm.Set_MixingAngle(0,1, 33.44 * Radian)
nsq_atm.Set_MixingAngle(0,2, 8.57 * Radian)
nsq_atm.Set_MixingAngle(1,2, T23BF)
nsq_atm.Set_SquareMassDifference(1,7.42e-5)
nsq_atm.Set_SquareMassDifference(2,DM31BF)
nsq_atm.Set_CPPhase(0,2,DCPBF)

nsq_atm.Set_rel_error(1.0e-4)
nsq_atm.Set_abs_error(1.0e-4)

nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
nsq_atm.EvolveState()



##################################### Event distribution    #####################################

Erec_min = 1
Erec_max = 10000
NErec = 40
erec = np.logspace(np.log10(Erec_min), np.log10(Erec_max), NErec+1, endpoint = True)
dlter = (np.log10(Erec_max) - np.log10(Erec_min))/NErec



Crec_min = -1
Crec_max = 1
Ncrec = 10
crec = np.linspace(Crec_min, Crec_max, Ncrec+1, endpoint = True)
dltcr = (Crec_max - Crec_min)/Ncrec


rate_weight = np.zeros_like(input_data["weight"])

mask_ErecMin = input_data["reco_energy"] > Erec_min

#Bin = np.zeros((2,len(input_data[mask_ErecMin])), dtype=int)
#Bin[0] = ((np.log10(input_data[mask_ErecMin]["reco_energy"]) - np.log10(Erec_min))/dlter).astype(int)
#Bin[1] = ((np.cos(input_data[mask_ErecMin]["reco_zenith"]) - Crec_min)/dltcr).astype(int)

#for i in range(len(input_data[mask_ErecMin])):
#    print(i, Bin[0][i])
#    if(input_data["reco_energy"][i] > Erec_min):
#        print(i, int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter), Bin[0][i])


mask_BinE = np.ones((NErec, len(input_data[mask_ErecMin]["reco_energy"])), dtype=bool)
for i in range(NErec) : 
    mask_BinE[i] = ((np.log10(input_data[mask_ErecMin]["reco_energy"]) - np.log10(Erec_min))/dlter).astype(int) == i 


mask_BinC = np.ones((Ncrec, len(input_data[mask_ErecMin]["reco_energy"])), dtype=bool)
for i in range(Ncrec) : 
    mask_BinC[i] = ((np.cos(input_data[mask_ErecMin]["reco_zenith"]) - Crec_min)/dltcr).astype(int) == i



Bin = np.zeros((2,len(input_data)), dtype=int)
for i in range(len(input_data)):
    if(input_data["reco_energy"][i] > Erec_min):
        Bin[0][i] = int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter)
        Bin[1][i] = int((np.cos(input_data["reco_zenith"][i]) - Crec_min)/dltcr)

        
EventsBF = np.zeros((2, NErec, Ncrec))
WEBF = np.zeros((2, len(rate_weight)))
WEBARBF = np.zeros((2, len(rate_weight)))
WMBF = np.zeros((2, len(rate_weight)))
WMBARBF = np.zeros((2, len(rate_weight)))


###################### Nue component BF  ############################################
nsq_atm.Set_initial_state(AtmIninue,nsq.Basis.flavor)
nsq_atm.EvolveState()
for i in range(len(rate_weight)):
#for i in range(len(rate_weight[mask_ErecMin])):

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


    #if input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:
    #WEBF[neuint][i] += input_data[mask_ErecMin]["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
    if input_data["reco_energy"][i] > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:
        WEBF[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq


###################### Nue bar component  ############################################
nsq_atm.Set_initial_state(AtmIninueBar,nsq.Basis.flavor)
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
        WEBARBF[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq




###################### Nu mu component BF  ############################################
nsq_atm.Set_initial_state(AtmIninum,nsq.Basis.flavor)
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
        WMBF[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq


                        
###################### Nu mu Bar component BF  ############################################
nsq_atm.Set_initial_state(AtmIninumBar,nsq.Basis.flavor)
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
        WMBARBF[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq



for s in range(2):
    for i in range(len(rate_weight)):
        EventsBF[s][Bin[0][i]][Bin[1][i]] += WEBF[s][i] + WEBARBF[s][i] + WMBF[s][i] + WMBARBF[s][i]
        
#H, xedges = np.histogram(input_data["reco_energy"]), bins=(erec), weights=WCBF)H, xedges = np.histogram(input_data["reco_energy"]), bins=(erec), weights=WCBF)
#


Events = np.zeros((2, NDM31, NT23, NErec, Ncrec))


WE = np.zeros((2, len(rate_weight)))
WEBAR = np.zeros((2, len(rate_weight)))
WM = np.zeros((2, len(rate_weight)))
WMBAR = np.zeros((2, len(rate_weight)))

for t in range(BT23, BT23+1):
    for m in range(BM31, BM31+1):
        for d in range(BDCP, BDCP+1):
#        for d in range(0, NDCP+1):

            ######################### Mixing parameters    ######################################
            nsq_atm.Set_MixingAngle(0,1, 33.44 * Radian)
            nsq_atm.Set_MixingAngle(0,2, 8.57 * Radian)
            nsq_atm.Set_MixingAngle(1,2, T23[t])
            nsq_atm.Set_SquareMassDifference(1,7.42e-5)
            nsq_atm.Set_SquareMassDifference(2,DM31[m])
            nsq_atm.Set_CPPhase(0,2,DCP[d])
            nsq_atm.Set_rel_error(1.0e-4)
            nsq_atm.Set_abs_error(1.0e-4)

            ################### Nue component     #################################
        
            nsq_atm.Set_initial_state(AtmIninue,nsq.Basis.flavor)
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
                    WE[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq


            ###################### Nue bar component  ############################################

            nsq_atm.Set_initial_state(AtmIninueBar,nsq.Basis.flavor)
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
                    WEBAR[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq




            ###################### Nu mu component  ############################################

            nsq_atm.Set_initial_state(AtmIninum,nsq.Basis.flavor)
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
                    WM[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq


                        
            ###################### Nu mu component  ############################################

            nsq_atm.Set_initial_state(AtmIninumBar,nsq.Basis.flavor)
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
                    WMBAR[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq



        


Xi = np.zeros((NDM31,NT23))


#############################              Systematics     ################################

#Normalzation
Norm_BF = 1.
SigNorm = 0.4


#Neutrino/anti-neutrino ratio
Bar_BF = 1.
SigBar = 0.02


#Flavor ratio
Flv_BF = 1.
SigFlv = 0.05


#Energy dependence of the flux
Gam_BF = 0.
SigGam = 0.2
E0Gam = 10.


#Horizontal vs Up or Down
Cos_BF = 0.
SigCos = 0.2


cosZen = np.cos(input_data["true_zenith"])
TanCos = np.tanh(cosZen) ** 2


#EventsE = np.zeros((2, NErec, Ncrec))
#EventsEBAR = np.zeros((2, NErec, Ncrec))
#EventsM = np.zeros((2, NErec, Ncrec))
#EventsMBAR = np.zeros((2, NErec, Ncrec))


#for i in range(len(rate_weight)):
#    if input_data["pid"][i] == 0 :
#        neuint = 0
#    elif input_data["pid"][i] == 1 :
#        neuint = 1
#
#    EventsE[neuint][Bin[0][i]][Bin[1][i]] += WE[neuint][i]
#    EventsEBAR[neuint][Bin[0][i]][Bin[1][i]] += WEBAR[neuint][i] 
#    EventsM[neuint][Bin[0][i]][Bin[1][i]] += WM[neuint][i]
#    EventsMBAR[neuint][Bin[0][i]][Bin[1][i]] += WMBAR[neuint][i]



mask_cUP = np.cos(input_data["true_zenith"]) < 0 
mask_cDOWN = np.cos(input_data["true_zenith"]) >= 0
#zenith = np.ones(len(input_data["true_energy"]))
#zenith[mask_cUP] = zenith[mask_cUP] - TanCos[mask_cUP] 
#print(len(zenith), zenith)


#exit(0)

def funXi(SYS):
    xxi = 0

    EventsE = np.zeros((2, NErec, Ncrec))
    EventsEBAR = np.zeros((2, NErec, Ncrec))
    EventsM = np.zeros((2, NErec, Ncrec))
    EventsMBAR = np.zeros((2, NErec, Ncrec))

    tilt = (input_data["true_energy"]/E0Gam) ** SYS[3]

    zenith = np.ones(len(input_data["true_energy"]))

    zenith[mask_cUP] = zenith[mask_cUP] - SYS[4] * TanCos[mask_cUP]
    zenith[mask_cDOWN] = zenith[mask_cDOWN] - SYS[5] * TanCos[mask_cDOWN]

    WEnew = np.zeros((2, len(rate_weight)))
    WEBARnew = np.zeros((2, len(rate_weight)))
    WMnew = np.zeros((2, len(rate_weight)))
    WMBARnew = np.zeros((2, len(rate_weight)))

    WEnew[0] = WE[0] * tilt * zenith
    WEnew[1] = WE[1] * tilt * zenith
    
    WEBARnew[0] = WEBAR[0] * tilt * zenith
    WEBARnew[1] = WEBAR[1] * tilt * zenith
    
    WMnew[0] = WM[0] * tilt * zenith
    WMnew[1] = WM[1] * tilt * zenith
    
    WMBARnew[0] = WMBAR[0] * tilt * zenith
    WMBARnew[1] = WMBAR[1] * tilt * zenith


    EventsE[0], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WEnew[0])
    EventsE[1], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WEnew[1])

    EventsEBAR[0], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WEBARnew[0])
    EventsEBAR[1], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WEBARnew[1])

    EventsM[0], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WMnew[0])
    EventsM[1], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WMnew[1])

    EventsMBAR[0], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WMBARnew[0])
    EventsMBAR[1], xedges, yedges = np.histogram2d(input_data["reco_energy"], np.cos(input_data["reco_zenith"]), bins=(erec,crec), weights=WMBARnew[1])


    #WEnew[0] = WE[0] * tilt
    #WEnew[1] = WE[1] * tilt
    #
    #WEBARnew[0] = WEBAR[0] * tilt
    #WEBARnew[1] = WEBAR[1] * tilt
    #
    #WMnew[0] = WM[0] * tilt
    #WMnew[1] = WM[1] * tilt
    #
    #WMBARnew[0] = WMBAR[0] * tilt
    #WMBARnew[1] = WMBAR[1] * tilt



    
    #for i in range(len(rate_weight)):
    #    if input_data["pid"][i] == 0 :
    #        neuint = 0
    #    elif input_data["pid"][i] == 1 :
    #        neuint = 1
    #    
    #    EventsE[neuint][Bin[0][i]][Bin[1][i]] += WEnew[neuint][i]
    #    EventsEBAR[neuint][Bin[0][i]][Bin[1][i]] += WEBARnew[neuint][i]
    #    EventsM[neuint][Bin[0][i]][Bin[1][i]] += WMnew[neuint][i]
    #    EventsMBAR[neuint][Bin[0][i]][Bin[1][i]] += WMBARnew[neuint][i]

        
    for s in range(2):
        for i in range(NErec) : 
            for j in range(Ncrec) : 

                if EventsBF[s][i][j] > 0 :
                    xxi += (SYS[0] * (SYS[2] * EventsE[s][i][j] + SYS[2] * SYS[1] * EventsEBAR[s][i][j] + EventsM[s][i][j] + SYS[1] * EventsMBAR[s][i][j]) - EventsBF[s][i][j]) * (SYS[0] * (SYS[2] * EventsE[s][i][j] + SYS[2] * SYS[1] * EventsEBAR[s][i][j] + EventsM[s][i][j] + SYS[1] * EventsMBAR[s][i][j]) - EventsBF[s][i][j]) / EventsBF[s][i][j]

    xxi += (SYS[0] - Norm_BF) * (SYS[0] - Norm_BF)/(SigNorm * SigNorm)    #normalization
    xxi += (SYS[1] - Bar_BF) * (SYS[1] - Bar_BF)/(SigBar * SigBar)        #neutrino/anti-neutrino
    xxi += (SYS[2] - Flv_BF) * (SYS[2] - Flv_BF)/(SigFlv * SigFlv)        #Flavor ratio
    xxi += (SYS[3] - Gam_BF) * (SYS[3] - Gam_BF)/(SigGam * SigGam)        #Tilt
    xxi += (SYS[4] - Cos_BF) * (SYS[4] - Cos_BF)/(SigCos * SigCos)        #Direction
    xxi += (SYS[5] - Cos_BF) * (SYS[5] - Cos_BF)/(SigCos * SigCos)        #Direction
                                               

    return xxi



MinXi = 0
name = "Xi_T23_{}_DM31_{}_DCP_{}_Sysv2.dat".format(BT23,BM31,BDCP)
my_file = open(name, 'w')
#for d in range(NT13) :



X0 = np.zeros((6))
X0[0] = 1.05      #Normalization
X0[1] = 1.01      #Neutrino vs anti-neutrino
X0[2] = 1.01      #Electron vs muon
X0[3] = 0.05      #Energy Tilt
X0[4] = 0.05      #Direction Up
X0[5] = 0.05      #Direction Down

bnds = ((0., 2.), (0., 2.), (0., 2.), (-1., 1), (-1., 1), (-1., 1))

    
res = minimize(funXi, X0, method='L-BFGS-B', bounds=bnds, options={'disp' : True, 'gtol' : 1e-2})

print(round(T23[BT23], 3), round(DM31[BM31], 5), round(DCP[d], 5), MinXi)
my_file.write("{},  {},  {},  {},  {},  {},  {},  {},  {},  {}\n".format(T23[BT23], DM31[BM31], DCP[BDCP], funXi(res.x), res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5]))

 



