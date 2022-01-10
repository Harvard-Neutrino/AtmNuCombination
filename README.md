# AtmNuCombination

For the time being the code is contained under the src folder and structured in Analysis and Simulation folders.

- src/Analysis: Currently performs standard 3-flavor neutrino oscillation analyses for the implemented experiments (IceCube Upgrade and Super-Kamiokande) and their combination. The code also includes some systematic uncertainties associated to the neutrino source and each experiment. 
An analysis is specified to the code via a xml file which contains all the information about the experiments, the neutrino sources and the systematic errors. An example of this can be found in src/Analysis/xmlAnalysis/AnalysisTemplate.xml.
In order to run some examples and get familiar with the code, the src/Analysis/run_example.sh contains a few quick examples. The main program for running the analysis is src/Analysis/runAnalysis.py and is run as follows:
python3 runAnalysis.py <analysis_xml_file> <output_file>
The results of the analysis are saved in <output_file> as a text file of columns. Further, these output files can be plotted using src/Analysis/PlotGlobalSens.py. Some examples are shown in src/Analysis/plot_example.sh.
python3 PlotGlobalSens.py <experiment> <output_file> [output_file2]
In order to run src/Analysis/PlotGlobalSens.py, one should specify the experiment: SK, IC or IC+SK.
The option of providing a second analysis output file is only for the cases when one wants to add the output of two separate analyses but with the same parameter grid.




The folder utils gathers tools for converting ROOT files produced by GENIE to HDF5 files, and read and merge any HDF5 files.