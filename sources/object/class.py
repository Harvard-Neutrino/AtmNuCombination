from enum import Enum 
import numpy as np

# First define some enum classes
class Topology(Enum):
    cascade: 0
    track: 1

class NeuType(Enum):
    AntiNeutrino: -1
    Neutrin0: 1

class Flavor(Enum):
    e: 12
    mu: 14
    tau: 16

# Now define the class of a single MCevent
class MCEvent:
    def __init__(self, neutype, flavor, topology, t_energy, \
                t_zenith, r_energy, r_zenith, W_MC, flux = 0):
        self.neutype = NeuType(neutype)
        self.flavor = Flavor(flavor)
        self.topology = Topology(topology)
        self.true = {"true_energy": t_energy, "true_zenith": t_zenith}
        self.reco = {"reco_energy": r_energy, "reco_zenith": r_zenith}
        self.w_mc = W_MC
        self.flux = flux
    
    def set_flux(self, flux):
        self.flux = flux
    
    def view_event(self):
        print("Event Details:")
        # add event details here

# Now define the class of a MC simulation
# Each MC simulation class contains 8 arrays of events w/ different neutuype, flavor, top
# Each array is then binnable
class simulation:
    def __init__(self, input_file):
        self.input_file = input_file
        # initialize 12+1 arrays
        self.events = np.array([])
        self.nu_e_cascade = np.array([])
        self.nu_e_track = np.array([])
        self.nub_e_cascade = np.array([])
        self.nub_e_track = np.array([])
        self.nu_mu_cascade = np.array([])
        self.nu_mu_track = np.array([])
        self.nub_mu_cascade = np.array([])
        self.nub_mu_track = np.array([])

        # Now put the events into the arrays
        for i in range(len(input_file["true_energy"])):
            flavor = np.absolute(input_file["pdg"])
            neutype = input_data["pdg"] / flavor
            new_event = MCEvent(neutype, flavor, input_file["pid"][i], input_file["true_energy"][i], \
                                input_file["true_zenith"][i], input_file["reco_energy"][i], \
                                input_file["reco_zenith"][i], input_file["weight"][i])
            if new_event.flavor.value == 12:
                if new_event.neutype.value == 1:
                    if new_event.topology.value == 0:
                        np.append(self.nu_e_cascade, np.array(new_event), axis = 0)
                    elif new_event.topology.value == 1:
                        np.append(self.nu_e_track, np.array(new_event), axis = 0)
                elif new_event.neutype.value == -1:
                    if new_event.topology.value == 0:
                        np.append(self.nub_e_cascade, np.array(new_event), axis = 0)
                    elif new_event.topology.value == 1:
                        np.append(self.nub_e_track, np.array(new_event), axis = 0)
            elif new_event.flavor.value == 14:
                if new_event.neutype.value == 1:
                    if new_event.topology.value == 0:
                        np.append(self.nu_e_cascade, np.array(new_event), axis = 0)
                    elif new_event.topology.value == 1:
                        np.append(self.nu_e_track, np.array(new_event), axis = 0)
                elif new_event.neutype.value == -1:
                    if new_event.topology.value == 0:
                        np.append(self.nub_e_cascade, np.array(new_event), axis = 0)
                    elif new_event.topology.value == 1:
                        np.append(self.nub_e_track, np.array(new_event), axis = 0)
