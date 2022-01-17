import numpy as np
import classdef as cl 
from params import *
import sys
import os
import math

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

def num_id():
    return len(t23l) * len(m31l)

def id_to_2did(idx):
    t23id = math.floor(idx / len(m31l))
    m31id = idx % len(m31l)
    return t23id, m31id

def id_to_param(idx):
    t23id, m31id = id_to_2did(idx)
    t23val = t23l[t23id]
    m31val = m31l[m31id]
    return t23val, m31val 

def local():
    totidx = len(t23l) * len(m31l)
    print("index should be ", totidx)
    for i in range(totidx):
        file_name = "{}.txt".format(i)
        complete_name = os.path.join(dir_name, file_name)
        file1 = open(complete_name, "a")
        file1.close()

def read_output(dir_name = dir_name):
    res = np.ndarray(shape = (len(m31l), len(t23l)))
    # print(m31l)
    # print(res.shape)
    for filename in os.listdir(dir_name):
        if filename.endswith(".txt"):
            print(filename)
            with open(os.path.join(dir_name, filename), "r") as file1:
                try:
                    temp = [line.rstrip('\n') for line in file1]
                    # print(temp)
                    # print(temp)
                    # print(int(float(temp[0])))
                    # i = 0
                    # Lines = file1.readLines()
                    # for line in Lines:
                    #     temp[i] = float(line)
                    thidx, dmidx = id_to_2did(int(float(temp[0])))
                    res[dmidx][thidx] = float(temp[1])
                except:
                    # print("excepted")
                    index = int(filename.split(".")[0])
                    # print(index)
                    thidx, dmidx = id_to_2did(int(float(index)))
                    res[dmidx][thidx] = 0
                    continue
    print(res)
    return res 