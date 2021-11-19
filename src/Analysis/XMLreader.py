import xml.etree.ElementTree as ET

# def parseXML(xmlfile='test.xml'):
def parseXML(xmlfile='AnalysisTemplate.xml'):
	
	# create element tree object
	tree = ET.parse(xmlfile)
  
	# get root element
	root = tree.getroot()

	print('Neutrino sources considered:')
	# for item in root.find('NeutrinoSource'):
	for source in root:
		sourceON = int(source.find('status').text)
		if sourceON:
			print(' + ',source.attrib['name'])
			for syst in source:
				systON = int(source.find('status').text)
				stag = syst.tag
				if str(stag) == 'systematic' and systON:
					nuisance = float(syst.find('nuisance').text)
					isFree = not bool(int(syst.find('tofit').text))
					if isFree:
						print('   - ',syst.attrib['name'],' with nominal nuisance parameter ', nuisance, ' and let free in the analysis.')
					else:
						print('   - ',syst.attrib['name'],' with nuisance parameter ', nuisance, ' fixed in the analysis.')



				# systON = int(syst.find('status').text)
				# print(systON)
				# if systON:
				# 	print('    ', syst.attrib['name'])






	# print(root.get('Analysis/NeutrinoSource/Atmospheric/name'))
	# print([elem.tag for elem in root.iter()])
	# for item in root:
	# 	print(item.tag)
	# 	if item.tag == 'NeutrinoSource':
	# 		for r in item:
	# 			for rr in r:
	# 				print(rr.tag)
	# 				if str(rr.tag) == 'name':
	# 					name = str(rr.attrib)
	# 				elif rr.tag == 'status':
	# 					status == str(rr.attrib)
	# 			print('Status of ', name, ' is ', status)
  
	# # create empty list for news items
	# newsitems = []
  
	# # iterate news items
	# for item in root.findall('./channel/item'):
  
		# # empty news dictionary
		# news = {}
  
		# # iterate child elements of item
		# for child in item:
		  
			# # special checking for namespace object content:media
			# if child.tag == '{http://search.yahoo.com/mrss/}content':
				# news['media'] = child.attrib['url']
			# else:
				# news[child.tag] = child.text.encode('utf8')
			  
			# # append news dictionary to news items list
			# newsitems.append(news)
			  
	# # return news items list
	# return newsitems

parseXML()
