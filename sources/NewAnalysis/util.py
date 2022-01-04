import numpy as np
import classdef as cl 

def bundle_fluxes(cth_nodes, energy_nodes, theta23, dm31, dcp):
	bf_e_flux = cl.Flux(cth_nodes, energy_nodes)
	bf_e_flux.set_initial_flux(cl.Flavor.e, cl.NeuType.Neutrino)
	bf_e_flux.propagate_flux(theta23, dm31, dcp)

	bf_mu_flux = cl.Flux(cth_nodes, energy_nodes)
	bf_mu_flux.set_initial_flux(cl.Flavor.mu, cl.NeuType.Neutrino)
	bf_mu_flux.propagate_flux(theta23, dm31, dcp)

	bf_ebar_flux = cl.Flux(cth_nodes, energy_nodes)
	bf_ebar_flux.set_initial_flux(cl.Flavor.e, cl.NeuType.AntiNeutrino)
	bf_ebar_flux.propagate_flux(theta23, dm31, dcp)

	bf_mubar_flux = cl.Flux(cth_nodes, energy_nodes)
	bf_mubar_flux.set_initial_flux(cl.Flavor.mu, cl.NeuType.AntiNeutrino)
	bf_mubar_flux.propagate_flux(theta23, dm31, dcp)

	bf_fluxes = np.array([[bf_e_flux, bf_ebar_flux],[bf_mu_flux, bf_mubar_flux]])

	return bf_fluxes

def get_all_weights(analysis, pointtype):
	analysis.get_weights(cl.Flavor.e, cl.NeuType.Neutrino, pointtype)
	analysis.get_weights(cl.Flavor.mu, cl.NeuType.Neutrino, pointtype)
	analysis.get_weights(cl.Flavor.e, cl.NeuType.AntiNeutrino, pointtype)
	analysis.get_weights(cl.Flavor.mu, cl.NeuType.AntiNeutrino, pointtype)

	return