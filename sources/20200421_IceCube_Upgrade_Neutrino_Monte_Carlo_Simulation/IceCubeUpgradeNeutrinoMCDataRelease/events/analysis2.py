import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq
import nuflux

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})

# Define path to file (you may need to change this to match your system)
input_file = "neutrino_mc.csv"

# Load the file using pandas
input_data = pd.read_csv(input_file)


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
nue_nc_mask = nue_mask & nc_mask
numu_nc_mask = numu_mask & nc_mask
nutau_nc_mask = nutau_mask & nc_mask 

# Define some energy bins (used throughout this notebook)
energy_bins_fine = np.logspace(0., 2., num=21)
energy_bins_course = np.logspace(0., 2., num=11)


units = nsq.Const()

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


nsq_atm.Set_initial_state(AtmInitialFlux,nsq.Basis.flavor)
nsq_atm.Set_ProgressBar(True) # progress bar will be printed on terminal
nsq_atm.EvolveState()

lifetime = 5.0*365*24*60*60
meter_to_cm_sq = 1e4
rate_weight = np.zeros_like(input_data["weight"])
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
        
    if input_data["true_energy"][i]*units.GeV < E_min or input_data["true_energy"][i]*units.GeV > E_max:
        rate_weight[i] = 0
        continue
    rate_weight[i] = input_data["weight"][i]*nsq_atm.EvalFlavor(neuflavor,
                                                                np.cos(input_data["true_zenith"][i]),
                                                                input_data["true_energy"][i]*\
                                                                units.GeV,neutype)*lifetime*meter_to_cm_sq

'''
This next part compares different methods of calculating the flux
'''
#################################################################################
# The example flux                                                              #
# Function to generate a toy atmopsheric neutrino flux                          #   
def atmo_nu_flux(true_energy) :                                                 #
    numu_flux = 5e2 * np.power(true_energy, -3.5)                               #
    nue_flux = numu_flux / 2.                                                   #
    return nue_flux, numu_flux                                                  #
#################################################################################


#############################################################################################
# The single-energy mode calculated flux                                                    #
def get_rate_single(nu_energy, nu_cos_zenith, pdg, weight):                                 #
	# Set the propagation                                                               #
	nu = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)                                     #
	nu.Set_E(nu_energy*units.GeV)  	                                                    #
	nu.Set_Body(nsq.EarthAtm())                                                         #
	nu.Set_Track(nsq.EarthAtm.Track(np.arccos(nu_cos_zenith)))                          #
	# Grab the flux                                                                     #
	nueflux = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue                    #
	#print(nueflux)                                                                     #
	#nuebflux = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue               #
	numuflux = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu                 #
	#numubflux = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar        #
	nutauflux = flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau              #
	#nutaubflux = flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar     #
	# Set the initial Flux                                                              #
	nu.Set_initial_state(np.array([nueflux, numuflux, nutauflux]), nsq.Basis.flavor)    #
	#print(nu.EvalFlavor(0))                                                            #
	nu.EvolveState()                                                                    #
	# Get the oscillated flux                                                           #
	oscillatedNuE = nu.EvalFlavor(0)                                                    #
	oscillatedNuMu = nu.EvalFlavor(1)                                                   #
	oscillatedNuTau = nu.EvalFlavor(2)                                                  #
	#print(oscillatedNuE, oscillatedNuMu, oscillatedNuTau)                              #
	# now get the rate                                                                  #
	if pdg == 12:                                                                       #
		return weight * oscillatedNuE                                               #
	if pdg == 14:                                                                       #
		return weight * oscillatedNuMu                                              #
	if pdg == 16:                                                                       #
		return weight * oscillatedNuTau                                             #
	else:                                                                               #
		print("Not a Neutrino Event Selected")                                      #
		return oscillatedNuE, oscillatedNuMu, oscillatedNuTau                       #
#############################################################################################


#########################################################################
# multiple-energy mode calculated flux                                  #
def get_rate_multiple(i):                                               #
    return rate_weight[i]                                               #
#########################################################################



def plot_rate_true_energy(rate_weight):
    # First multiply by lifetime and conversion
    #rate_weight[i] *= lifetime*meter_to_cm_sq
    input_data["rate_weight"] = rate_weight
    # Note that converting to mHz for the rate, as this is a more suitable unit for the IceCube Upgrade 
    fig, ax = plt.subplots(figsize=(7,6))
    fig.suptitle("True Energy Rated Weight")
    ax.hist(input_data["true_energy"][nue_cc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][nue_cc_mask], \
          label=r"$\nu_{e,CC}$", color="blue", histtype="step")
    ax.hist(input_data["true_energy"][numu_cc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][numu_cc_mask], \
          label=r"$\nu_{\mu,CC}$", color="red", histtype="step")
    ax.hist(input_data["true_energy"][nutau_cc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][nutau_cc_mask], \
          label=r"$\nu_{\tau,CC}$", color="green", histtype="step")
    ax.hist(input_data["true_energy"][nc_mask], bins=energy_bins_fine, \
          weights=(1e3)*input_data["rate_weight"][nc_mask], \
          label=r"$\nu_{NC}$", color="grey", histtype="step")
    
    print("Total neutrino rate = %0.3g mHz" % (np.nansum(rate_weight) * 1e3) )
    
    ax.set_xscale("log")
    ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
    ax.set_xlim(10, 100)
    ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
                     useOffset=None, useLocale=None, useMathText=None)
    ax.set_ylabel("Rate [mHz]")
    ax.grid(True)
    _ = ax.legend()
    fig.savefig("rated_weight_distribution_true_energy(mHz).png")

# plot_rate_true_energy(rate_weight)
	
	
	
def plot_rate_reconstructed_energy(rate_weight):
    # First multiply by lifetime and conversion
    #rate_weight[i] *= lifetime*meter_to_cm_sq
    input_data["rate_weight"] = rate_weight
    # Note that converting to mHz for the rate, as this is a more suitable unit for the IceCube Upgrade 
    fig, ax = plt.subplots(figsize=(7,6))
    fig.suptitle("Reconstructed Energy Rated Weight")
    ax.hist(input_data["reco_energy"][nue_cc_mask], bins=energy_bins_fine, \
          weights=(1e-4)*input_data["rate_weight"][nue_cc_mask], \
          label=r"$\nu_{e,CC}$", color="blue", histtype="step")
    ax.hist(input_data["reco_energy"][numu_cc_mask], bins=energy_bins_fine, \
          weights=(1e-4)*input_data["rate_weight"][numu_cc_mask], \
          label=r"$\nu_{\mu,CC}$", color="red", histtype="step")
    ax.hist(input_data["reco_energy"][nutau_cc_mask], bins=energy_bins_fine, \
          weights=(1e-4)*input_data["rate_weight"][nutau_cc_mask], \
          label=r"$\nu_{\tau,CC}$", color="green", histtype="step")
    ax.hist(input_data["reco_energy"][nc_mask], bins=energy_bins_fine, \
          weights=(1e-4)*input_data["rate_weight"][nc_mask], \
          label=r"$\nu_{NC}$", color="grey", histtype="step")
    
    print("Total neutrino rate = %0.3g mHz" % (np.nansum(rate_weight) * 1e3) )
    
    ax.set_xscale("log")
    ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
    ax.set_xlim(10, 100)
    ax.ticklabel_format(axis='y', style='sci', scilimits=None,\
                     useOffset=None, useLocale=None, useMathText=None)
    ax.set_ylabel("Rate [mHz]")
    ax.grid(True)
    ax.legend()
    fig.savefig("rated_weight_distribution_reco_energy(mHz).png")

	
# plot_rate_reconstructed_energy(rate_weight)

def plot_rate_comparison_true_energy(rate_weight):
	input_data["rate_weight"] = rate_weight
	
	# Get cumulative rates
	nueCC = 1e3 * np.nansum(input_data["rate_weight"][nue_cc_mask])
	nueNC = 1e3 * np.nansum(input_data["rate_weight"][nue_nc_mask])
	numuCC = 1e3 * np.nansum(input_data["rate_weight"][numu_cc_mask])
	numuNC = 1e3 * np.nansum(input_data["rate_weight"][numu_nc_mask])
	nutauCC = 1e3 * np.nansum(input_data["rate_weight"][nutau_cc_mask])
	nutauNC = 1e3 * np.nansum(input_data["rate_weight"][nutau_nc_mask])
	
	currents = ['CC', 'NC']
	nue = np.array([nueCC, nueNC])
	numu = np.array([numuCC, numuNC])
	nutau = np.array([nutauCC, nutauNC])
	ind = [x for x, _ in enumerate(currents)]
	
	plt.subplots(figsize=(7,6))
	plt.bar(ind, nue, width=0.8, label=r"$\nu_{e}$", color='skyblue', bottom=numu+nutau)
	plt.bar(ind, numu, width=0.8, label=r"$\nu_{\mu}$", color='lightcoral', bottom=nutau)
	plt.bar(ind, nutau, width=0.8, label=r"$\nu_{\tau}$", color='seagreen')

	plt.xticks(ind, currents)
	plt.ylabel("Weighted Event Rates (mHz)")
	plt.xlabel("Interaction Types")
	plt.legend(loc="upper right")
	plt.title("Weighted Event Rates Make Up", y = 1.08)
	plt.ticklabel_format(axis='y', style='sci', scilimits=None,\
                     useOffset=None, useLocale=None, useMathText=None)
	plt.savefig("Weighted_Event_Rates_Make_Up(True_Energy).png")

plot_rate_comparison_true_energy(rate_weight)
	














