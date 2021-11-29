import numpy as np
import math
from math import asin, sqrt
from SimReader import Reader, OscPar
import multiprocessing
import sys
from XMLreader import parseXML

# Setup experiment details
analysis_xml_file=str(sys.argv[1])
outfile=str(sys.argv[2])

an = parseXML(analysis_xml_file)
an.readSources()
an.readExperiments()
an.readSystematics()
an.readPhysics()
an.readOscPar()

mc = Reader(an.experiments[0],an.mcFiles[0])
mc.Binning()
mc.SetOutFile(outfile)

# Get unoscillated atm. fluxes
mc.InitialFlux()

# Setup oscillation parameters grid and best fit value
mc.BFOscillatorxml(an.neutrinos,**an.OscParametersBest)

manager = multiprocessing.Manager()
return_dict = manager.dict()

processes = []
j = 0

with open(outfile,'w') as f:
	f.write('SqT12 SqT13 SqT23 DM21 DM31 DCP Ordering X2\n')
# Osc. in parameters space (multiprocessing)
for i,t12 in enumerate(an.OscParametersGrid['Sin2Theta12']):
	for j,t13 in enumerate(an.OscParametersGrid['Sin2Theta13']):
		for k,t23 in enumerate(an.OscParametersGrid['Sin2Theta23']):
			for l,dm21 in enumerate(an.OscParametersGrid['Dm221']):
				for m,dm31 in enumerate(an.OscParametersGrid['Dm231']):
					for d,cp in enumerate(an.OscParametersGrid['dCP']):
						for o,hi in enumerate(an.OscParametersGrid['Ordering']):
							p = multiprocessing.Process(target=mc.Chi2Calculator,args=[an.neutrinos,t12, t13, t23, dm21, dm31, cp, hi])
							if __name__ == "__main__":
								processes.append(p)
								p.start()
								j = j + 1
								print(f'{t12} {t13} {t23} {dm21} {dm31} {cp} process started')
								if j%20==0:
									for i,p in enumerate(processes):
										p.join()
									processes = []
