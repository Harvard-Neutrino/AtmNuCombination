
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq
import nuflux

from sterileprob import prob		

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'patch.linewidth': 3})

# Define path to file (you may need to change this to match your system)
input_file = "neutrino_mc.csv"

# Load the file using pandas
input_data = pd.read_csv(input_file)

print("input data loaded")

# Define units
units = nsq.Const()

# Define some energy bins (used throughout this notebook)
energy_bins_fine = np.logspace(0., 2., num=21)
energy_bins_course = np.logspace(0., 2., num=11)

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

print("Defined bins and masks")

# Calc effective area
def eff_area():
	effective_area_hist_ecc, bin_edges_ecc = np.histogram(input_data["true_energy"][nue_cc_mask],\
	                                                  weights=input_data["weight"][nue_cc_mask],\
	                                                  bins=energy_bins_fine )
	effective_area_hist_ecc /= 4. * np.pi # Normalise by solid angle (using the fully sky)
	effective_area_hist_ecc /= np.diff(bin_edges_ecc) # Bin widths

	effective_area_hist_mucc, bin_edges_mucc = np.histogram(input_data["true_energy"][numu_cc_mask],\
	                                                  weights=input_data["weight"][numu_cc_mask],\
	                                                  bins=energy_bins_fine )
	effective_area_hist_mucc /= 4. * np.pi # Normalise by solid angle (using the fully sky)
	effective_area_hist_mucc /= np.diff(bin_edges_mucc) # Bin widths

	effective_area_hist_taucc, bin_edges_taucc = np.histogram(input_data["true_energy"][nutau_cc_mask],\
	                                                  weights=input_data["weight"][nutau_cc_mask],\
	                                                  bins=energy_bins_fine )
	effective_area_hist_taucc /= 4. * np.pi # Normalise by solid angle (using the fully sky)
	effective_area_hist_taucc /= np.diff(bin_edges_taucc) # Bin widths

	effective_area_hist_nc, bin_edges_nc = np.histogram(input_data["true_energy"][nc_mask],\
	                                                  weights=input_data["weight"][nc_mask],\
	                                                  bins=energy_bins_fine )
	effective_area_hist_nc /= 4. * np.pi # Normalise by solid angle (using the fully sky)
	effective_area_hist_nc /= np.diff(bin_edges_nc) # Bin widths

	# Plot
	fig, ax = plt.subplots(figsize=(7,5))
	ax.step( bin_edges_ecc, effective_area_hist_ecc.tolist()+[effective_area_hist_ecc[-1]],\
	        where="post", color="blue", label=r"$\nu_{e,CC}$" )
	ax.step( bin_edges_mucc, effective_area_hist_mucc.tolist()+[effective_area_hist_mucc[-1]],\
	        where="post", color="red", label=r"$\nu_{\mu,CC}$" )
	ax.step( bin_edges_taucc, effective_area_hist_taucc.tolist()+[effective_area_hist_taucc[-1]],\
	        where="post", color="green", label=r"$\nu_{\tau,CC}$" )
	ax.step( bin_edges_nc, effective_area_hist_nc.tolist()+[effective_area_hist_nc[-1]],\
	        where="post", color="gray", label= "All " r"$\nu_{NC}$" )
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
	ax.set_ylabel(r"$A_{\rm{eff}}$ [$\rm{m^2}$]")
	ax.grid(True)
	_ = ax.legend()

	#fig.savefig("EffectiveArea.png")

neutrino_flavors = 3

flux = nuflux.makeFlux("honda2006")

# First lets check out what the flux looks like
def plot_flux():
	fig, ax = plt.subplots(figsize = (7,6))
	energy_grid, coszen_grid = np.meshgrid(energy_nodes, cth_nodes, indexing="ij")
	print(coszen_grid)
	flux_grid = flux.getFlux(nuflux.NuE, energy_grid/units.GeV, np.arccos(coszen_grid))
	cmesh = ax.pcolormesh(energy_grid, coszen_grid, flux_grid, vmin=0., vmax=1000.)
	fig.colorbar(cmesh, ax=ax, label=r"$P(\nu_e)$"" flux")
	ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
	ax.set_xscale("log")
	ax.set_ylabel(r"$\cos(\theta_{\rm{zenith,true}})$")
	#fig.savefig("Nu_E_Fluxx.png")

# This gets the oscillated flux of one event
def get_rate(nu_energy, nu_cos_zenith, pdg, weight):
	# Set the propagation
	nu = nsq.nuSQUIDS(3, nsq.NeutrinoType.neutrino)
	nu.Set_E(nu_energy*units.GeV)
	nu.Set_Body(nsq.EarthAtm())
	nu.Set_Track(nsq.EarthAtm.Track(np.arccos(nu_cos_zenith)))
	# Grab the flux
	nueflux = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
	#print(nueflux)
	#nuebflux = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue
	numuflux = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
	#numubflux = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
	nutauflux = flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau
	#nutaubflux = flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar
	# Set the initial Flux
	nu.Set_initial_state(np.array([nueflux, numuflux, nutauflux]), nsq.Basis.flavor)
	#print(nu.EvalFlavor(0))
	nu.EvolveState()
	# Get the oscillated flux
	oscillatedNuE = nu.EvalFlavor(0)
	oscillatedNuMu = nu.EvalFlavor(1)
	oscillatedNuTau = nu.EvalFlavor(2)
	#print(oscillatedNuE, oscillatedNuMu, oscillatedNuTau)
	# now get the rate
	if pdg == 12:
		return weight * oscillatedNuE
	if pdg == 14:
		return weight * oscillatedNuMu
	if pdg == 16:
		return weight * oscillatedNuTau
	else:
		print("Not a Neutrino Event Selected")
		return oscillatedNuE, oscillatedNuMu, oscillatedNuTau

#get_rate(50941380148/units.GeV, -0.48718, 12, 1)

def plot_rate():
	print("plotting rated weight distribution")
	# Plot histograms of event rates vs energy
	# First use get_rate function to obtain rate
	rate_weight = np.zeros_like(input_data["weight"])
	# for i in nue_mask:
	# 	rate_weight[i] = get_rate(input_data["true_energy"][i], \
	# 									input_data["true_zenith"][i], 12, \
	# 									input_data["weight"][i])
	# for i in numu_mask:
	# 	rate_weight[i] = get_rate(input_data["true_energy"][i], \
	# 									input_data["true_zenith"][i], 14, \
	# 									input_data["weight"][i])
	# for i in nutau_mask:
	# 	rate_weight[i] = get_rate(input_data["true_energy"][i], \
	# 									input_data["true_zenith"][i], 16, \
	# 									input_data["weight"][i])
	# input_data["rate_weight"] = rate_weight

	for i in range(len(rate_weight)):
		neutype = np.abs(input_data["pdg"][i])
		if neutype not in [12, 14, 16]:
			rate_weight[i] = 0
		else:
			rate_weight[i] = get_rate(input_data["true_energy"][i], \
										input_data["true_zenith"][i], \
										input_data["pdg"][i], \
										input_data["weight"][i])
	input_data["rate_weight"] = rate_weight

	# Note that converting to mHz for the rate, as this is a more suitable unit for the IceCube Upgrade 
	fig, ax = plt.subplots(figsize=(7,6))
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
	ax.set_xscale("log")
	ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
	ax.set_ylabel("Rate [mHz]")
	ax.grid(True)
	_ = ax.legend()
	fig.savefig("rated_weight_distribution_true_energy.png")


plot_rate()

def plot_rate_2(e, cth, pdg, weight):
	# first get the flux
	def getflux(nu_energy, nu_cos_zenith):
		nueflux = flux.getFlux(nuflux.NuE,nu_energy,nu_cos_zenith) # nue
		#nuebflux = flux.getFlux(nuflux.NuEBar,nu_energy,nu_cos_zenith) # nue
		numuflux = flux.getFlux(nuflux.NuMu,nu_energy,nu_cos_zenith) # numu
		#numubflux = flux.getFlux(nuflux.NuMuBar,nu_energy,nu_cos_zenith) # numu bar
		nutauflux = flux.getFlux(nuflux.NuTau,nu_energy,nu_cos_zenith) # nutau
		#nutaubflux = flux.getFlux(nuflux.NuTauBar,nu_energy,nu_cos_zenith) # nutau bar
		return nueflux, numuflux, nutauflux
	nueflux, numuflux, nutauflux, = getflux(e, cth)
	# then get the oscillation probabilities
	# Use 3+1 model of oscillation, mu to tau oscillation probability
	def prob(true_energy, true_zenith):
	    # find baseline length from true zenith
	    true_coszen = np.cos(true_zenith)
	    earth_radius = 6371.
	    production_height = 15. # Assuming neutrino produced 15 km above surface
	    detector_depth = 1. # Assuming detector depth of 1 km
	    baseline = -earth_radius*true_coszen +  \
	                np.sqrt( (earth_radius*true_coszen)**2 - earth_radius**2 + \
	                        (earth_radius+production_height+detector_depth)**2 )
    
	    # now compute the oscillation probability
	    Petau = prob(1, 3, "t2k", "IH", "t2k", manE = True, valE = true_energy, \
	              manL = True, valL = baseline)
	    Pemu = prob(1, 2, "t2k", "IH", "t2k", manE = True, valE = true_energy, \
	              manL = True, valL = baseline)
	    Pmutau = prob(2, 3, "t2k", "IH", "t2k", manE = True, valE = true_energy, \
	              manL = True, valL = baseline)
	    Pmue = prob(2, 1, "t2k", "IH", "t2k", manE = True, valE = true_energy, \
	              manL = True, valL = baseline)
	    Ptaue = prob(3, 1, "t2k", "IH", "t2k", manE = True, valE = true_energy, \
	              manL = True, valL = baseline)
	    Ptaumu = prob(3, 2, "t2k", "IH", "t2k", manE = True, valE = true_energy, \
	              manL = True, valL = baseline)
	    return Pemu, Petau, Pmue, Pmutau, Ptaue, Ptaumu

	# Get fluxes
	nue_flux, numu_flux, nutau_flux = getflux(input_data["true_energy"], input_data["true_zenith"])

	# Get oscillations
	Pemu, Petau, Pmue, Pmutau, Ptaue, Ptaumu = prob(e, cth)

	# Now compute weights representing event rates
	# This is flux * osc * weight
	rate_weight = np.zeros_like(input_data["weight"])
	rate_weight[nue_mask] = nue_flux[nue_mask] * input_data["weight"][nue_mask] # No oscillations in 2 flavor model
	rate_weight[numu_mask] = numu_flux[numu_mask] * (1. - Pmutau[numu_mask]) * input_data["weight"][numu_mask] # numu flux, less those that oscillation
	rate_weight[nutau_mask] = numu_flux[nutau_mask] * Pmutau[nutau_mask] * input_data["weight"][nutau_mask] # nutau purely come from oscillations in the numu flux
	input_data["rate_weight_new"] = rate_weight







