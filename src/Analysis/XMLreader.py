import xml.etree.ElementTree as ET
import sys
import numpy as np

class parseXML:
	def __init__(self, xmlfile='AnalysisTemplate.xml'):
		# create element tree object
		self.xmlfile = xmlfile
		self.tree = ET.parse(xmlfile)
		# get root element of XML file
		self.root = self.tree.getroot()
		self.disabledsyst = []

	def reader(self, item, atrib='name'):
		itemList = []
		for source in self.root.iter(item):
			if atrib=='name':
				if int(source.find('status').text):
					itemList.append(source.attrib[atrib])
				else:
					for elt in source.iter('systematic'):
						self.disabledsyst.append(elt.attrib['name'])
			else:
				itemList.append(source.attrib[atrib])
		return itemList

	def readSources(self):
		self.sources = self.reader('NeutrinoSource')
		print('------------------------------------')
		print('Neutrino sources considered:')
		for s in self.sources:
			print(' + ',s)
		print('====================================')

	def readExperiments(self):
		self.experiments = self.reader('NeutrinoExperiment')
		self.mcFiles = self.reader('simulation','filename')
		print('------------------------------------')
		print('Experiments considered:')
		for s,m in zip(self.experiments,self.mcFiles):
			print(' + ',s, ', at ',m)
		print('====================================')

	def readSystematics(self):
		self.systematics = self.reader('systematic')
		self.systematics = [x for x in self.systematics if x not in self.disabledsyst]
		print('------------------------------------')
		print('Systematic uncertainties considered:')
		for s in self.systematics:
			print(' + ',s)
		print('====================================')

	def readPhysics(self):
		self.physics = self.reader('NeutrinoPhysics')
		print('------------------------------------')
		print('Neutrino physics considered:')
		for s in self.physics:
			print(' + ',s)
		print('====================================')
	
	def readOscPar(self):
		keys = ['Sin2Theta12','Sin2Theta13','Sin2Theta23','Dm221','Dm231','dCP','Ordering']
		self.parameters = keys
		self.OscParametersGrid = {}
		self.OscParametersBest = {}
		for node in self.root.iter('NeutrinoPhysics'):
			for phys in self.physics:
				if node.attrib['name'] == phys:
					self.neutrinos = int(node.find('flavours').text)
					for key in keys:
						par = node.find('parameters/'+key)
						self.OscParametersGrid[key] = self.readOscParValues(par, key)
						self.OscParametersBest[key] = self.readOscParBestValue(par, key)

	def readOscParValues(self, pars, key):
		if key=='Ordering':
			if int(pars.find('normal').text):
				x = np.array(['normal'])
				if int(pars.find('inverted').text):
					x = np.append((x,'inverted'))
			else:
				if int(pars.find('inverted').text):
					x = np.array(['inverted'])
			return x


			return np.array([int(pars.find('normal').text),int(pars.find('inverted').text)])
		else:
			n = int(pars.find('points').text)
			mini = float(pars.find('min').text)
			maxi = float(pars.find('max').text)
			return np.linspace(mini, maxi, n, endpoint = True)

	def readOscParBestValue(self, pars, key):
		if key=='Ordering':
			return str(pars.find('best').text)
		else:
			return float(pars.find('best').text)


import itertools
if __name__ == '__main__':

	analysis = parseXML('AnalysisTemplate.xml')
	analysis.readSources()
	analysis.readExperiments()
	analysis.readSystematics()
	analysis.readPhysics()
	analysis.readOscPar()

	for element in itertools.product([analysis.OscParametersGrid[key] for key in analysis.parameters]):
		print(element)
	


