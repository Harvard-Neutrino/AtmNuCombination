# AtmNuCombination

For the time being the code is contained under the src folder and structured in Analysis and Simulation folders.

## Analysis
### src/Analysis
Currently performs standard 3-flavor neutrino oscillation analyses for the implemented experiments (IceCube Upgrade and Super-Kamiokande) and their combination. The code also includes some systematic uncertainties associated to the neutrino source and each experiment. 

### Dependencies
- numpy  
- nuSQuIDS  
- nuflux  
- h5py  
- matplotlib  
- pandas  
- scipy 

Except for nuSQuIDS, you can install them by doing: 
```
pip install -r requirements_analysis.txt
```
For nuSQuIDS, please follow the instructions at https://github.com/arguelles/nuSQuIDS/ .

### Running
An analysis is specified to the code via a xml file which contains all the information about the experiments, the neutrino sources and the systematic errors. An example of this can be found in src/Analysis/xmlAnalysis/AnalysisTemplate.xml.  

In order to run some examples and get familiar with the code, the src/Analysis/run_example.sh contains a few quick examples. The main program for running the analysis is src/Analysis/runAnalysis.py and is run as follows:  
```
  usage: runAnalysis.py [-h] [-p [POINT]] [-o [OUTFILE]] [--multi] [--cluster] [xml_file]  

  positional arguments:  
    xml_file              Input analysis file in xml format  
  
  optional arguments:  
    -h, --help            show this help message and exit  
    -p [POINT], --point [POINT]  
                          Specify analysis point to run. Only if 'cluster' option is enabled  
    -o [OUTFILE], --outfile [OUTFILE]  
                          Analysis output file  
    --multi               Option for running the analysis with multiprocessing (recommended locally)  
    --cluster             Option for submitting jobs to a cluster  
```

  
The results of the analysis are saved in <output_file> as a text file of columns (by deafult the output file is out.dat). If no running mode (--multi or --cluster) is specified it will run sequentially the list of points from the xml file.   

### Plotting
Further, these output files can be plotted using src/Analysis/PlotGlobalSens.py. Some examples are shown in src/Analysis/plot_example.sh.
```

python3 PlotGlobalSens.py <experiment> <output_file> [output_file2]

```

In order to run src/Analysis/plotting/PlotGlobalSens.py, one should specify the experiment: SK, IC or IC+SK.

The option of providing a second analysis output file is only for the cases when one wants to add the output of two separate analyses but with the same parameter grid.

NOTE: This needs to be improved and merged with Ivan's plotting.

## Simulation
### src/Simulations
 It contains the official MC simulation of IceCube Upgrade and preliminary effective/toy simulations for the ORCA and Super-Kamiokande experiments.

### Requirements
- nuSQuIDS
- nuSQUIDSTools
- nuflux
- particle
- matplotlib
- numpy
- pythia8
- h5py
- pandas
- scipy

#### Super-Kamiokande
Contains the relevant information extracted from various publications to closely match the reconstruction performance for SK. The code works applying the reconstruction on GENIE simulation files with format _gst_.

##### Running

```
Super-Kamiokande atmospheric neutrino simulation
usage: makeSimulation_HDF5.py [-h] [-o OUTFILENAME] [-v] [--sk] [--H] [--Gd] [in_hdf5filename]

positional arguments:
  in_hdf5filename       Input Genie HDF5 file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILENAME, --outfilename OUTFILENAME
                        Output file name in HDF5 format.
  -v, --verbo           Verbosity of simulation process.
  --SK                  Default SK simulation without neutron tagging.
  --H                   SKIV simulation with neutron tagging on hydrogen.
  --Gd                  SKVII simulation with neutron tagging on gadolinium.

```

#### IceCube Upgrage
Public release of the MC simulation at https://icecube.wisc.edu/data-releases/2020/04/icecube-upgrade-neutrino-monte-carlo-simulation/ .

#### ORCA
Tune IC-Up MC release to match the quoted resolution and efficiencies reported by the ORCA collaboration.

## Utils
The folder utils gathers tools for converting ROOT files produced by GENIE to HDF5 files, and read and merge any HDF5 files.

## Questions?
Please, send any doubts or suggestions to pablo.fernandez@dipc.org, ivan.j.martinez-soler@durham.ac.uk, carguelles@g.harvard.edu, miaochenjin@g.harvard.edu, santiagoginer@college.harvard.edu .
