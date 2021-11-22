import xml.etree.ElementTree as ET

class parseXML:
	def __init__(self, xmlfile='AnalysisTemplate.xml'):
		# create element tree object
		self.xmlfile = xmlfile
		self.tree = ET.parse(xmlfile)
		# get root element of XML file
		self.root = self.tree.getroot()
		self.disabledsyst = []

	def reader(self, item, atrib):
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
		self.sources = self.reader('NeutrinoSource','name')
		print('------------------------------------')
		print('Neutrino sources considered:')
		for s in self.sources:
			print(' + ',s)
		print('====================================')

	def readExperiments(self):
		self.experiments = self.reader('NeutrinoExperiment','name')
		self.MCFiles = self.reader('simulation','filename')
		print('------------------------------------')
		print('Experiments considered:')
		for s,m in zip(self.experiments,self.MCFiles):
			print(' + ',s, ', at ',m)
		print('====================================')

	def readSystematics(self):
		self.systematics = self.reader('systematic','name')
		self.systematics = [x for x in self.systematics if x not in self.disabledsyst]
		print('------------------------------------')
		print('Systematic uncertainties considered:')
		for s in self.systematics:
			print(' + ',s)
		print('====================================')


	def readPhysics(self):
		# self.sources = self.reader('NeutrinoParameters')
		pass
	


if __name__ == '__main__':

	analysis = parseXML('AnalysisTemplate.xml')
	analysis.readSources()
	analysis.readExperiments()
	analysis.readSystematics()
	# analysis.readPhysics()
