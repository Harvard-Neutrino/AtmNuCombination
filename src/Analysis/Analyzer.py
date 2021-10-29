import numpy as np
import pandas as pd
import h5py
#import nuSQUIDSpy as nsq
import nuSQuIDS as nsq
import nuflux
import sys
from scipy.optimize import minimize

units = nsq.Const()

Time = 3.0*365*24*60*60
meter_to_cm_sq = 1e4
Radian = np.pi / 180.


Bin = int(sys.argv[1])

BT23 = int(Bin/50)
BM31 = int(Bin - BT23 * 50)

experiment = str(sys.argv[2])


########################    READ DATA   #####################################

# Define path to file (you may need to change this to match your system)
if experiment=='SuperK':
    input_file = '../Simulation/SuperK/data/output/combined.hdf5'
    # input_file = 'combined.hdf5'
    SKNorm = 0.05411673990568942 # Bit of magic
    with h5py.File(input_file,'r') as hf:
        d_evis = np.array(hf['evis'][()])
        d_recocz = np.array(hf['recodirZ'][()])
        d_truecz = np.array(hf['dirnuZ'][()])
        d_azi = np.array(hf['azi'][()])
        d_mode = np.array(hf['mode'][()])
        d_ipnu = np.array(hf['ipnu'][()])
        d_pnu = np.array(hf['pnu'][()])
        d_weightReco = np.array(hf['weightReco'][()])
        d_weightSim = np.array(hf['weightSim'][()])
        d_itype = np.array(hf['itype'][()])
    condition1 = (d_itype<14) * (d_itype>-1) 
    condition2 = (d_itype<16) * (d_itype>13) * (d_evis>1)
    condition = (condition2 + condition1) * (d_evis<400)
    evis = d_evis[condition]
    recocz = d_recocz[condition]
    truecz = d_truecz[condition]
    azi = d_azi[condition]
    mode = d_mode[condition]
    ipnu = d_ipnu[condition]
    pnu = d_pnu[condition]
    weightReco = d_weightReco[condition]
    weightSim = d_weightSim[condition]
    itype = d_itype[condition]
else:
    input_file = "neutrino_mc.csv"
    # Load the file using pandas
    input_data = pd.read_csv(input_file)



#######################  USE DATA    #####################################

if experiment=='SuperK':
    # Define some reco energy bins for each SK sample (used throughout this notebook)
    sge_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
    sgm_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
    sgsrpi0ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 500])
    sgmrpi0ebins = np.array([0.1, 0.15, 0.25, 0.4, 0.63, 500])
    mge_ebins = np.array([1.3, 2.5, 5., 10., 500.])
    mgm_ebins = np.array([1.3, 3.0, 500.])
    mre_ebins = np.array([1.3, 2.5, 5.0, 500.])
    mrm_ebins = np.array([0.6, 1.3, 2.5, 5., 500.])
    mro_ebins = np.array([1.3, 2.5, 5.0, 10., 500.])
    pcs_ebins = np.array([0.1, 20., 1.0e5])
    pct_ebins = np.array([0.1, 10.0, 50., 1.0e3, 1.0e5])
    z10bins = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    z1bins = np.array([-1, 1.0])
    energy_bins = {0:sge_ebins, 1:sge_ebins, 2:sgsrpi0ebins, 3:sgm_ebins, 4:sgm_ebins, 5:sgm_ebins, 6:sgmrpi0ebins,
    7:mge_ebins, 8:mge_ebins, 9:mgm_ebins, 10:mre_ebins, 11:mre_ebins, 12:mrm_ebins, 13:mro_ebins, 14:pcs_ebins, 15:pct_ebins}
    cz_bins = {0:z10bins, 1:z1bins, 2:z1bins, 3:z10bins, 4:z10bins, 5:z1bins, 6:z1bins,
    7:z10bins, 8:z10bins, 9:z10bins, 10:z10bins, 11:z10bins, 12:z10bins, 13:z10bins, 14:z10bins, 15:z10bins}

    # Define masks to identify different neutrino flavors
    nue_mask = (np.abs(ipnu) == 12)
    numu_mask = (np.abs(ipnu) == 14)
    nutau_mask = (np.abs(ipnu) == 16)

    # Define masks to identify different flavor/interaction combinations.
    nc_mask = (np.abs(mode)>30)
    cc_mask = (np.abs(mode)<30)
    nue_cc_mask = nue_mask & cc_mask
    numu_cc_mask = numu_mask & cc_mask
    nutau_cc_mask = nutau_mask & cc_mask


else:
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

NT23 = 50
SqT23_min = 0.305
SqT23_max = 0.705
SqT23BF = np.sin(49.2*Radian) * np.sin(49.2*Radian)
SqT23 = np.linspace(SqT23_min, SqT23_max, NT23+1, endpoint = True)
T23 = np.arcsin(np.sqrt(SqT23))
T23BF = T23[33]



NDM31 = 50
DM31_max = 2.8e-3
DM31_min = 2.3e-3
DM31 = np.linspace(DM31_min, DM31_max, NDM31+1, endpoint = True)
DM31BF = DM31[22]


##################    Flux     ###########################

if experiment == 'SuperK':
    flux = nuflux.makeFlux('IPhonda2014_sk_solmin')
    E_min = 0.1
    E_max = 4.0e2
else:
    flux = nuflux.makeFlux('IPhonda2014_spl_solmin')
    E_min = 1.0
    E_max = 1.0e3


interactions = False

E_nodes = 100
energy_range = nsq.logspace(E_min,E_max,E_nodes)
energy_nodes = nsq.logspace(E_min*units.GeV,E_max*units.GeV,E_nodes)

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
for ic,nu_cos_zenith in enumerate(cth_nodes):
    for ie,nu_energy in enumerate(energy_range):
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

nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
nsq_atm.EvolveState()



##################################### Event distribution    #####################################

if experiment != 'SuperK':
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


if experiment=='SuperK':
    weight = weightSim*weightReco
    rate_weight = np.zeros_like(weight)
else:
    rate_weight = np.zeros_like(input_data["weight"])



if experiment=='SuperK':
    Bin = np.zeros((2,len(rate_weight)), dtype=int)
    # Asign a bin number to each event
    for i in range(len(rate_weight)):
        h, bins = np.histogram(evis[i],energy_bins[itype[i]])
        index = np.where(h==1)
        Bin[0][i] = int(index[0][0])
        h, bins = np.histogram(recocz[i],cz_bins[itype[i]])
        index = np.where(h==1)
        Bin[1][i] = int(index[0][0])
        # print(Bin[0][i])
        # print(Bin[1][i])
    EventsBF = np.zeros((16, 5, 10))
else:
    Bin = np.zeros((2,len(rate_weight)), dtype=int)
    # Asign a bin number to each event
    for i in range(len(rate_weight)):
        if input_data["reco_energy"][i] * units.GeV > Erec_min:
            Bin[0][i] = int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter)
            Bin[1][i] = int((np.cos(input_data["reco_zenith"][i]) - Crec_min)/dltcr)
    EventsBF = np.zeros((2, NErec, Ncrec))
    EvBFEne = np.zeros((2, NErec))
    EvBFCzen = np.zeros((2, Ncrec))


for i in range(len(rate_weight)):
    if experiment=='SuperK':
        if ipnu[i]>0:
            neutype=0
        else:
            neutype=1

        neuint = int(itype[i])

        if np.abs(ipnu[i])==12:
            neuflavor=0
        elif np.abs(ipnu[i])==14:
            neuflavor=1
        elif np.abs(ipnu[i])==16:
            neuflavor=2
        # EventsBF[neuint][Bin[0][i]][Bin[1][i]] += weight[i] * nsq_atm.EvalFlavor(neuflavor, truecz[i], pnu[i] * units.GeV, neutype) * SKNorm
        EventsBF[neuint][Bin[0][i]][Bin[1][i]] += weight[i] * 1 * SKNorm

    else:
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

            #binr = int( (np.log10(input_data["reco_energy"][i]) - np.log10(Erec_min))/dlter)
            #binc = int((np.cos(input_data["reco_zenith"][i]) - Crec_min)/dltcr)

            EventsBF[neuint][Bin[0][i]][Bin[1][i]] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
            
            EvBFEne[neuint][Bin[0][i]] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

            EvBFCzen[neuint][Bin[1][i]] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq



if experiment=='SuperK':
    Events = np.zeros((16, NDM31, NT23, 5, 10))
    EventsMass = np.zeros((16, NDM31, 5, 10))
    EventsMix = np.zeros((16, NT23, 5, 10))

    WE = np.zeros((16, len(rate_weight)))
    WEBAR = np.zeros((16, len(rate_weight)))
    WM = np.zeros((16, len(rate_weight)))
    WMBAR = np.zeros((16, len(rate_weight)))

else:
    Events = np.zeros((2, NDM31, NT23, NErec, Ncrec))
    #EventsE = np.zeros((2, NDM31, NT23, NErec, Ncrec))
    #EventsEBAR = np.zeros((2, NDM31, NT23, NErec, Ncrec))
    #EventsM = np.zeros((2, NDM31, NT23, NErec, Ncrec))
    #EventsMBAR = np.zeros((2, NDM31, NT23, NErec, Ncrec))
    EventsMass = np.zeros((2, NDM31, NErec, Ncrec))
    EventsMix = np.zeros((2, NT23, NErec, Ncrec))


    WE = np.zeros((2, len(rate_weight)))
    WEBAR = np.zeros((2, len(rate_weight)))
    WM = np.zeros((2, len(rate_weight)))
    WMBAR = np.zeros((2, len(rate_weight)))


for t in range(BT23, BT23+1):
    for m in range(BM31, BM31+1):


        ######################### Mixing parameters    ######################################
        nsq_atm.Set_MixingAngle(0,1, 33.44 * Radian)
        nsq_atm.Set_MixingAngle(0,2, 8.57 * Radian)
        nsq_atm.Set_MixingAngle(1,2, T23[t])
        nsq_atm.Set_SquareMassDifference(1,7.42e-5)
        nsq_atm.Set_SquareMassDifference(2,DM31[m])

        if experiment=='SuperK':
                    ################### Nue component     #################################
                        
            nsq_atm.Set_initial_state(AtmIninue,nsq.Basis.flavor)
            nsq_atm.EvolveState()
            
            for i in range(len(rate_weight)):

                if ipnu[i] > 0 :
                    neutype = 0
                else:
                    neutype = 1
                    
                neuint = int(itype[i])

                if np.abs(ipnu[i]) == 12:
                    neuflavor = 0
                elif np.abs(ipnu[i]) == 14:
                    neuflavor = 1
                elif np.abs(ipnu[i]) == 16:
                    neuflavor = 2

                # WE[neuint][i] += weight[i] * nsq_atm.EvalFlavor(neuflavor, truecz[i], pnu[i] * units.GeV, neutype) * SKNorm
                WE[neuint][i] += weight[i] * 1 * SKNorm


            ###################### Nue bar component  ############################################
            
            nsq_atm.Set_initial_state(AtmIninueBar,nsq.Basis.flavor)
            nsq_atm.EvolveState()
            
            for i in range(len(rate_weight)):

                if ipnu[i] > 0 :
                    neutype = 0
                else:
                    neutype = 1

                neuint = int(itype[i]) 

                if np.abs(ipnu[i]) == 12:
                    neuflavor = 0
                elif np.abs(ipnu[i]) == 14:
                    neuflavor = 1
                elif np.abs(ipnu[i]) == 16:
                    neuflavor = 2

                # WEBAR[neuint][i] += weight[i] * nsq_atm.EvalFlavor(neuflavor, truecz[i], pnu[i] * units.GeV, neutype) * SKNorm
                WEBAR[neuint][i] += weight[i] * 1 * SKNorm


            ###################### Nu mu component  ############################################
            
            nsq_atm.Set_initial_state(AtmIninum,nsq.Basis.flavor)
            nsq_atm.EvolveState()
            
            for i in range(len(rate_weight)):

                if ipnu[i] > 0 :
                    neutype = 0
                else:
                    neutype = 1

                neuint = int(itype[i]) 

                if np.abs(ipnu[i]) == 12:
                    neuflavor = 0
                elif np.abs(ipnu[i]) == 14:
                    neuflavor = 1
                elif np.abs(ipnu[i]) == 16:
                    neuflavor = 2

                # WM[neuint][i] += weight[i] * nsq_atm.EvalFlavor(neuflavor, truecz[i], pnu[i] * units.GeV, neutype) * SKNorm
                WM[neuint][i] += weight[i] * 1 * SKNorm


            ###################### Nu mu component  ############################################
            
            nsq_atm.Set_initial_state(AtmIninumBar,nsq.Basis.flavor)
            nsq_atm.EvolveState()
            
            for i in range(len(rate_weight)):

                if ipnu[i] > 0 :
                    neutype = 0
                else:
                    neutype = 1

                neuint = int(itype[i]) 

                if np.abs(ipnu[i]) == 12:
                    neuflavor = 0
                elif np.abs(ipnu[i]) == 14:
                    neuflavor = 1
                elif np.abs(ipnu[i]) == 16:
                    neuflavor = 2

                # WMBAR[neuint][i] += weight[i] * nsq_atm.EvalFlavor(neuflavor, truecz[i], pnu[i] * units.GeV, neutype) * SKNorm
                WMBAR[neuint][i] += weight[i] * 1 * SKNorm

                    

        else:
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


                if input_data["reco_energy"][i] * units.GeV > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:


                    #EventsE[neuint][m][t][Bin[0][i]][Bin[1][i]] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

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


                if input_data["reco_energy"][i] * units.GeV > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:

                    #EventsEBAR[neuint][m][t][Bin[0][i]][Bin[1][i]] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
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


                if input_data["reco_energy"][i] * units.GeV > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:

                    #EventsM[neuint][m][t][Bin[0][i]][Bin[1][i]] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
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


                if input_data["reco_energy"][i] * units.GeV > Erec_min  and input_data["true_energy"][i] * units.GeV > E_min and input_data["true_energy"][i] * units.GeV < E_max:

                    #EventsMBAR[neuint][m][t][Bin[0][i]][Bin[1][i]] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq
                    WMBAR[neuint][i] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq


                    

                #Events[neuint][m][t][binr][binc] += input_data["weight"][i] * nsq_atm.EvalFlavor(neuflavor, np.cos(input_data["true_zenith"][i]), input_data["true_energy"][i] * units.GeV, neutype) * Time * meter_to_cm_sq

        



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

TanCos = np.zeros((len(rate_weight)))
for i in range(len(rate_weight)):
    if experiment=='SuperK':
        cosZen = truecz[i]
    else:
        cosZen = np.cos(input_data["true_zenith"][i])
    TanCos[i] = np.tanh(cosZen) ** 2

MinXi = 1e10

if experiment=='SuperK':
    EventsE = np.zeros((16, 5, 10))
    EventsEBAR = np.zeros((16, 5, 10))
    EventsM = np.zeros((16, 5, 10))
    EventsMBAR = np.zeros((16, 5, 10))
else:   
    EventsE = np.zeros((2, NErec, Ncrec))
    EventsEBAR = np.zeros((2, NErec, Ncrec))
    EventsM = np.zeros((2, NErec, Ncrec))
    EventsMBAR = np.zeros((2, NErec, Ncrec))


def funXi(SYS):
    xxi = 0

    if experiment=='SuperK':
        sge_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
        sgm_ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 1.33])
        sgsrpi0ebins = np.array([0.1, 0.25, 0.4, 0.63, 1.0, 500])
        sgmrpi0ebins = np.array([0.1, 0.15, 0.25, 0.4, 0.63, 500])
        mge_ebins = np.array([1.3, 2.5, 5., 10., 500.])
        mgm_ebins = np.array([1.3, 3.0, 500.])
        mre_ebins = np.array([1.3, 2.5, 5.0, 500.])
        mrm_ebins = np.array([0.6, 1.3, 2.5, 5., 500.])
        mro_ebins = np.array([1.3, 2.5, 5.0, 10., 500.])
        pcs_ebins = np.array([0.1, 20., 1.0e5])
        pct_ebins = np.array([0.1, 10.0, 50., 1.0e3, 1.0e5])
        z10bins = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        z1bins = np.array([-1, 1.0])
        energy_bins = {0:sge_ebins, 1:sge_ebins, 2:sgsrpi0ebins, 3:sgm_ebins, 4:sgm_ebins, 5:sgm_ebins, 6:sgmrpi0ebins,
        7:mge_ebins, 8:mge_ebins, 9:mgm_ebins, 10:mre_ebins, 11:mre_ebins, 12:mrm_ebins, 13:mro_ebins, 14:pcs_ebins, 15:pct_ebins}
        cz_bins = {0:z10bins, 1:z1bins, 2:z1bins, 3:z10bins, 4:z10bins, 5:z1bins, 6:z1bins,
        7:z10bins, 8:z10bins, 9:z10bins, 10:z10bins, 11:z10bins, 12:z10bins, 13:z10bins, 14:z10bins, 15:z10bins}
        for s in range(16):
            for i in range(5) : 
                for j in range(10) : 
                    EventsE[s][i][j] = 0.
                    EventsEBAR[s][i][j] = 0.
                    EventsM[s][i][j] = 0.
                    EventsMBAR[s][i][j] = 0.


        for i in range(len(rate_weight)):
            if ipnu[i] == 0 :
                neuint = 0
            elif ipnu[i] == 1 :
                neuint = 1
            
            tilt = (pnu[i] / E0Gam) ** SYS[3]


            if cosZen <= 0 :
                tcos =  1 - SYS[4] * TanCos[i]
            else :
                tcos =  1 - SYS[5] * TanCos[i]

            EventsE[neuint][Bin[0][i]][Bin[1][i]] += WE[neuint][i] * tilt * tcos
            EventsEBAR[neuint][Bin[0][i]][Bin[1][i]] += WEBAR[neuint][i] * tilt * tcos
            EventsM[neuint][Bin[0][i]][Bin[1][i]] += WM[neuint][i] * tilt * tcos
            EventsMBAR[neuint][Bin[0][i]][Bin[1][i]] += WMBAR[neuint][i] * tilt * tcos

            
        for s in range(16):
            for i in range(energy_bins[s].size-1) : 
                for j in range(cz_bins[s].size-1) : 
                    
                    if EventsBF[s][i][j] > 0 :
                        xxi += (SYS[0] * (SYS[2] * EventsE[s][i][j] + SYS[2] * SYS[1] * EventsEBAR[s][i][j] + EventsM[s][i][j] + SYS[1] * EventsMBAR[s][i][j]) - EventsBF[s][i][j]) * (SYS[0] * (SYS[2] * EventsE[s][i][j] + SYS[2] * SYS[1] * EventsEBAR[s][i][j] + EventsM[s][i][j] + SYS[1] * EventsMBAR[s][i][j]) - EventsBF[s][i][j]) / EventsBF[s][i][j]
                        
                                 
        xxi += (SYS[0] - Norm_BF) * (SYS[0] - Norm_BF)/(SigNorm * SigNorm)    #normalization
        xxi += (SYS[1] - Bar_BF) * (SYS[1] - Bar_BF)/(SigBar * SigBar)        #neutrino/anti-neutrino
        xxi += (SYS[2] - Flv_BF) * (SYS[2] - Flv_BF)/(SigFlv * SigFlv)        #Flavor ratio
        xxi += (SYS[3] - Gam_BF) * (SYS[3] - Gam_BF)/(SigGam * SigGam)        #Tilt
        xxi += (SYS[4] - Cos_BF) * (SYS[4] - Cos_BF)/(SigCos * SigCos)        #Direction
        xxi += (SYS[5] - Cos_BF) * (SYS[5] - Cos_BF)/(SigCos * SigCos)        #Direction

    else:
        for s in range(2):
            for i in range(NErec) : 
                for j in range(Ncrec) : 
                    EventsE[s][i][j] = 0.
                    EventsEBAR[s][i][j] = 0.
                    EventsM[s][i][j] = 0.
                    EventsMBAR[s][i][j] = 0.


        for i in range(len(rate_weight)):
            if input_data["pid"][i] == 0 :
                neuint = 0
            elif input_data["pid"][i] == 1 :
                neuint = 1
            
            tilt = (input_data["true_energy"][i] / E0Gam) ** SYS[3]


            if cosZen <= 0 :
                tcos =  1 - SYS[4] * TanCos[i]
            else :
                tcos =  1 - SYS[5] * TanCos[i]

            EventsE[neuint][Bin[0][i]][Bin[1][i]] += WE[neuint][i] * tilt * tcos
            EventsEBAR[neuint][Bin[0][i]][Bin[1][i]] += WEBAR[neuint][i] * tilt * tcos
            EventsM[neuint][Bin[0][i]][Bin[1][i]] += WM[neuint][i] * tilt * tcos
            EventsMBAR[neuint][Bin[0][i]][Bin[1][i]] += WMBAR[neuint][i] * tilt * tcos

            
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




#for t in range(BT23, BT23+1):
#    for m in range(BM31, BM31+1):
#    #for m in range(NDM31):
#
#        MinXi = 1e10
#
#        for n in range(NNorm):
#
#                xxi = 0
#
#                for s in range(2):
#                        for i in range(NErec) : 
#                                for j in range(Ncrec) : 
#                
#                                        if EventsBF[s][i][j] > 0 :
#                                                xxi += (Norm[n] * (EventsE[s][m][t][i][j] + EventsEBAR[s][m][t][i][j] + EventsM[s][m][t][i][j] + EventsMBAR[s][m][t][i][j]) - EventsBF[s][i][j]) * (Norm[n] * (EventsE[s][m][t][i][j] + EventsEBAR[s][m][t][i][j] + EventsM[s][m][t][i][j] + EventsMBAR[s][m][t][i][j]) - EventsBF[s][i][j]) / EventsBF[s][i][j]
#                                                #xxi += (Norm[n] * Events[s][m][t][i][j] - EventsBF[s][i][j]) * (Norm[n] * Events[s][m][t][i][j] - EventsBF[s][i][j]) / EventsBF[s][i][j]
#                             
#                xxi += (Norm[n] - 1) * (Norm[n] -1 )/(SigNorm * SigNorm)
#
#                if xxi < MinXi :
#                        MinXi = xxi
#
#        Xi[m][t] = MinXi


X0 = np.zeros((6))

X0[0] = 1.05      #Normalization
X0[1] = 1.01      #Neutrino vs anti-neutrino
X0[2] = 1.01      #Electron vs muon
X0[3] = 0.05      #Energy Tilt
X0[4] = 0.05      #Direction Up
X0[5] = 0.05      #Direction Down

bnds = ((0., 2.), (0., 2.), (0., 2.), (-1., 1), (-1., 1), (-1., 1))

res = minimize(funXi, X0, method='L-BFGS-B', bounds=bnds, options={'disp' : True})

print(round(T23[BT23], 3), round(DM31[BM31], 5))

name = "Xi_T23_{}_DM31_{}_Sys.dat".format(BT23,BM31)
np.savetxt(name, np.c_[round(T23[BT23], 3), round(DM31[BM31], 5), funXi(res.x), res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5]], delimiter=' ')
#np.savetxt(name, np.c_[round(T23[BT23], 3), round(DM31[BM31], 10), funXi(res.x), res.x], delimiter=' ')
#np.savetxt(name, np.c_[round(T23[BT23], 3), round(DM31[BM31], 2), Xi[BM31][BT23], res.x], delimiter=' ')
 
for t in range(BT23, BT23+1):
    for m in range(BM31, BM31+1):
    #for m in range(NDM31):
        print(T23[t], DM31[m], funXi(res.x), res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5])



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