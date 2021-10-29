import xml.etree.ElementTree as ET

def parseXML(xmlfile='../AnalysisTemplate.xml'):
	
	# create element tree object
	tree = ET.parse(xmlfile)
  
	# get root element
	root = tree.getroot()

	status = 'dummy'
	name = 'dummy'
	[elem.tag for elem in root.iter()]
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
