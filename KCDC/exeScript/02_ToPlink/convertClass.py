'''
Created on 2016. 1. 27.

@author: HMY
'''

from fileInOut import fileInOut
fileIO = fileInOut()

class convertClass(object):
	print "CLASS: convertClass"
	
	
	def covertToMAP(self, callName, callIn, annoIn, mapOut):	
		print "CLASS: Running coverToMAP()..."
		
		annoData = fileIO.fileRead(annoIn)
		callData = fileIO.fileRead(callIn)
		
		annoDic = {}
		for i in range(1, len(annoData)):
			annoArray = annoData[i].split()
			annoDic.setdefault(annoArray[0], annoArray[3] + ' ' + annoArray[0] + ' 0 ' + annoArray[4])
		
		if callName[0:17] == "AxiomGT1.calls.0_":
			iSTART = 1
		else:
			iSTART = 0
		
		mapWrite = open(mapOut, 'w')
		for i in range(iSTART, len(callData)):
			callArray = callData[i].split('\t')
			callID = callArray[0]
			
			if callID in annoDic:
				mapWrite.write(annoDic[callID] + '\n')
		mapWrite.close()
	
	
	def convertToPED(self, callName, callIn, sampleIn, pedOut):
		print "CLASS: Running converToPED()..."
		
		callData = fileIO.fileRead(callIn)
		sampleData = fileIO.fileRead(sampleIn)
		
		if callName[0:17] == "AxiomGT1.calls.0_":
			iSTART = 1
		else:
			iSTART = 0
		
		for i in range(int(iSTART), len(callData)):
			callArray = callData[i].split()
			for j in range(1, len(callArray)):
				callGeno = int(callArray[j])
				if callGeno == 0:
					annoGeno = "A A"
				elif callGeno == 1:
					annoGeno = "A B"
				elif callGeno == 2:
					annoGeno = "B B"
				elif callGeno == -1:
					annoGeno = "0 0"
				else:
					print(annoGeno)
				#print len(callArray), i, j
				sampleData[j] = sampleData[j] + ' ' + annoGeno
		
		pedWrite = open(pedOut, 'w')
		for i in range(1, len(sampleData)):
			pedWrite.write(sampleData[i] + '\n')
		pedWrite.close()
			
