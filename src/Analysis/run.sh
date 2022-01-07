#!/bin/bash

time python3 runAnalysis.py xmlAnalysis/SK_T23_DM31_stats.xml SK/SK_T23_DM31_stats.dat
time python3 runAnalysis.py xmlAnalysis/IC_T23_DM31_stats.xml IC/IC_T23_DM31_stats.dat
time python3 runAnalysis.py xmlAnalysis/IC+SK_T23_DM31_stats.xml IC+SK/IC+SK_T23_DM31_stats.dat
