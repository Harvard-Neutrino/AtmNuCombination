import sys
from SimReader import Reader
from xmlReader import parseXML
from Sensitivity import sensitivity
import argparse

# Read arguments
############################
parser = argparse.ArgumentParser()
parser.add_argument("xml_file", type=str, nargs='?', default='xmlAnalysis/AnalysisTemplate.xml', help='Input analysis file in xml format. Just used for getting the true values.')
parser.add_argument('-o', '--outfile', nargs='?', type=str, default='out.dat', help='Analysis output file.')
args = parser.parse_args()

# Setup analysis files
############################
analysis_xml_file = args.xml_file
outfile = args.outfile

# Setup analysis from xml file
############################
an = parseXML(analysis_xml_file)
an.readSources()
an.readExperiments()
an.readDetectors()
an.readPhysics()
an.readOscPar()
an.CheckSystematics()

# Setup all experiments
############################
mcList = {}
for s in an.sources:
	for i,(exp,fil,t) in enumerate(zip(an.experiments,an.mcFiles,an.Exposure)):
		# print(s,exp,t,fil)
		mcList[exp] = Reader(s,exp,t,fil)
		mcList[exp].Binning()
		# Get unoscillated atm. fluxes
		mcList[exp].InitialFlux()
		# Set best fit value oscillations
		mcList[exp].BFOscillator(an.neutrinos,**an.OscParametersBest)

print('=============================================================')
print('==================== Starting analysis ======================')
print('=============================================================')

""" We're only interested in comparing both orderings keeping the rest of the 
parameters the same (may not be true in more in-depth analyses). """

OscParametersBest_alternateMO = an.OscParametersBest
if an.OscParametersBest['Ordering'] == 'normal':
	OscParametersBest_alternateMO['Ordering'] = 'inverted'
else:
	OscParametersBest_alternateMO['Ordering'] = 'normal'

binned_sensitivity(OscParametersBest_alternateMO[:-1], OscParametersBest_alternateMO[-1], an, mcList, outfile)


