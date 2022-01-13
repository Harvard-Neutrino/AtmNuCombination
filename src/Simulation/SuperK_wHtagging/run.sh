#!/bin/bash

# Original SK
for i in data/input/atm_emu/gntp*hdf5; do python3 makeSimulation_HDF5.py $i; done
for i in data/input/atm_tau/gntp*hdf5; do python3 makeSimulation_HDF5.py $i; done

# SK with H-neutron tagging
for i in data/input/atm_emu/gntp*hdf5; do python3 makeSimulation_HDF5.py $i --H; done
for i in data/input/atm_tau/gntp*hdf5; do python3 makeSimulation_HDF5.py $i --H; done

# SK with Gd-neutron tagging
for i in data/input/atm_emu/gntp*hdf5; do python3 makeSimulation_HDF5.py $i --Gd; done
for i in data/input/atm_tau/gntp*hdf5; do python3 makeSimulation_HDF5.py $i --Gd; done

# Merge all simulated files
python3 ../../../utils/HDF5merger.py data/output/SK data/output/SK/combined.hdf
python3 ../../../utils/HDF5merger.py data/output/SK_Htag data/output/SK_Htag/combined.hdf
python3 ../../../utils/HDF5merger.py data/output/SK_Gdtag data/output/SK_Gdtag/combined.hdf

# Re-weighting after combination
python3 ammendWeight data/output/SK/combined.hdf
python3 ammendWeight data/output/SK_Htag/combined.hdf
python3 ammendWeight data/output/SK_Gdtag/combined.hdf