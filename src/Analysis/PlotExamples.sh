#!/bin/bash

# SK contours for NO
python3 GlobalSens.py SK SK/X2_stat_25_20_NO_NO.txt

# SK contours for both orderings
python3 GlobalSens.py SK SK/X2_stat_25_20_NO_NO.txt SK/X2_stat_25_20_NO_IO.txt

# ICUp contours for NO
python3 GlobalSens.py IC IC/Chi2_ST23_DM31_DCP_NoSys.dat

# SK and ICUp combination contours for NO
python3 GlobalSens.py IC+SK IC/Chi2_ST23_DM31_DCP_NoSys.dat SK/Chi2_ST23_DM31_DCP_NoSys.dat