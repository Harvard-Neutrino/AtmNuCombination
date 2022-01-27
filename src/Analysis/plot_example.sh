#!/bin/bash

# Plots the sensitvity plots for the analyses done in run_example.sh.
python3 PlotGlobalSens.py SK SK/SK_T23_DM31_stats.dat
python3 PlotGlobalSens.py IC IC/IC_T23_DM31_stats.dat
python3 PlotGlobalSens.py IC+SK IC+SK/IC+SK_T23_DM31_stats.dat

# The latter should be equal to 
# python3 GlobalSens.py IC+SK SK/SK_T23_DM31_stats.dat IC/IC_T23_DM31_stats.dat


