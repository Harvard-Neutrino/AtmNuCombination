#!/bin/bash

# Computes the standard 3-flavor oscillation analysis for SK and IC atm. neutrinos and their combination assuming no systematics.
python3 runAnalysis.py xmlAnalysis/SK_T23_DM31_stats.xml SK/SK_T23_DM31_stats.dat
python3 runAnalysis.py xmlAnalysis/IC_T23_DM31_stats.xml IC/IC_T23_DM31_stats.dat
python3 runAnalysis.py xmlAnalysis/IC+SK_T23_DM31_stats.xml IC+SK/IC+SK_T23_DM31_stats.dat
