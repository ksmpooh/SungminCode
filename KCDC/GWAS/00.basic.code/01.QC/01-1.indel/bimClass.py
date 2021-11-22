'''
Created on 2017. 11. 6.

@author: myhwang
'''

from fileInOut import fileInOut
fileIO = fileInOut()

class bimClass(object):
	print "CLASS: bimClass"
	
	
	def convertToBim(self, annoIn, bimIn, bimOut):
		print "Running convertToBim()..."
		
		annoData = fileIO.fileRead(annoIn)
		bimData = fileIO.fileRead(bimIn)
		
		annoDic = {}
		for i in range(1, len(annoData)):
			annoArray = annoData[i].split('\t')
			annoID = annoArray[0]
			data = annoArray[7] + ' ' + annoArray[8]
			annoDic.setdefault(annoID, data)
		
		bimWrite = open(bimOut, 'w')
		for i in bimData:
			bimArray = i.split()
			bimID = bimArray[1]
			bimREF = bimArray[4]
			bimALT = bimArray[5]
			
			if bimID in annoDic:
				annoArray = annoDic[bimID].split()
				Aallele = annoArray[0]
				Ballele = annoArray[1]
				
				if bimREF == "A" and bimALT == "A":
					bimWrite.write(' '.join(bimArray[:4]) + ' ' + Aallele + ' ' + Aallele + '\n')
				elif bimREF == "A" and bimALT == "B":
					bimWrite.write(' '.join(bimArray[:4]) + ' ' + Aallele + ' ' + Ballele + '\n')
				elif bimREF == "B" and bimALT == "B":
					bimWrite.write(' '.join(bimArray[:4]) + ' ' + Ballele + ' ' + Ballele + '\n')
				elif bimREF == "B" and bimALT == "A":
					bimWrite.write(' '.join(bimArray[:4]) + ' ' + Ballele + ' ' + Aallele + '\n')
				elif bimREF == "A" and bimALT == "0":
					bimWrite.write(' '.join(bimArray[:4]) + ' ' + Aallele + ' 0\n')
				elif bimREF == "B" and bimALT == "0":
					bimWrite.write(' '.join(bimArray[:4]) + ' ' + Ballele + ' 0\n')
				elif bimREF == "0" and bimALT == "A":
					bimWrite.write(' '.join(bimArray[:4]) + ' 0 ' + Aallele + '\n')
				elif bimREF == "0" and bimALT == "B":
					bimWrite.write(' '.join(bimArray[:4]) + ' 0 ' + Ballele + '\n')
				elif bimREF == "0" and bimALT == "0":
					bimWrite.write(' '.join(bimArray[:4]) + ' 0 0 \n')
		bimWrite.close()
