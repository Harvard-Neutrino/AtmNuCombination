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


#Define some plotting parameters for future use
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

nsq_atm = nsq.nuSQUIDSAtm(cth_nodes,energy_nodes,neutrino_flavors,nsq.NeutrinoType.both,interactions)

flux = nuflux.makeFlux("honda2006")


# First lets check out what the flux looks like
fig, ax = plt.subplots(figsize = (7,6))
energy_grid, coszen_grid = np.meshgrid(energy_nodes, cth_nodes, indexing="ij")
flux_grid = flux.getFlux(nuflux.NuE, energy_grid/units.GeV, np.arccos(coszen_grid))
cmesh = ax.pcolormesh(energy_grid, coszen_grid, flux_grid, vmin=0., vmax=1000.)
fig.colorbar(cmesh, ax=ax, label=r"$P(\nu_e)$"" flux")
ax.set_xlabel(r"$E_{\nu,\rm{true}}$ [GeV]")
ax.set_xscale("log")
ax.set_ylabel(r"$\cos(\theta_{\rm{zenith,true}})$")
fig.savefig("Nu_E_Fluxx.png")

'''
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

AtmInitialFluxNuE = np.zeros((len(cth_nodes), len(energy_nodes))) 
print(AtmInitialFluxNuE.shape)
for ic in range(0, len(cth_nodes)):
    for ie in range(0, len(energy_nodes)):
        AtmInitialFluxNuE = AtmInitialFlux[ic][ie][0][0]

print("HAAAAAAAAAAAAAAAAAAAAAAAA")
print(AtmInitialFluxNuE.shape)

plt.imshow(AtmInitialFluxNuE, cmap = "hot", interpolation = "nearest")
plt.savefig("ExampleFlux.png")
'''
