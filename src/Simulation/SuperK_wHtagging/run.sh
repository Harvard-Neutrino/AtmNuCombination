#!/bin/bash

for i in data/input/atm_emu/gntp*hdf5; do python3 makeSimulation_HDF5.py $i; done
for i in data/input/atm_tau/gntp*hdf5; do python3 makeSimulation_HDF5.py $i; done
