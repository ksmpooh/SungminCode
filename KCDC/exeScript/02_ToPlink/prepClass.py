'''
Created on 2016. 1. 27.

@author: HMY
'''

from fileInOut import fileInOut
fileIO = fileInOut()

class prepClass(object):
	print "CLASS: prepClass"
	
	
	def posRange(self, plinkDir, mergeOut, posOut, posSTART, posEND, posUnit):
		print "CLASS: Running posRange()..."
		
		cnt = 0
		mergeWrite = open(mergeOut , 'w')
		for i in range(posSTART, posEND, posUnit):
			if cnt == 0:
				start = str(i) + ','
				end = str(posUnit) + ','
				cnt = cnt + 1
			elif posEND <= (i + posUnit):
				start = start + str(i)
				end = end + "END"
				mergeWrite.write(plinkDir + "KNIH.RAW." + str(i) + "_END.ped " + 
								plinkDir + "KNIH.RAW." + str(i) + "_END.map\n")
			else:
				start = start + str(i) + ','
				end = end + str(i + posUnit) + ','
				mergeWrite.write(plinkDir + "KNIH.RAW." + str(i) + '_' + str(i + posUnit) + ".ped " + 
								plinkDir + "KNIH.RAW." + str(i) + '_' + str(i + posUnit) + ".map\n")
		
		posWrite = open(posOut, 'w')
		posWrite.write(start + '\n')
		posWrite.write(end + '\n')
		posWrite.close()
		mergeWrite.close()
		
		return start, end
	
	
	def sampleID(self, callIn, sexIn, sampleOut):
		print "CLASS: Running sampeID()..."
		
		callData = fileIO.fileRead(callIn)
		sexData = fileIO.fileRead(sexIn)
		
		sexDic = {}
		for i in sexData:
			sexArray = i.split()
			#sexDic.setdefault(sexArray[0], sexArray[1])
			sexDic.setdefault(sexArray[0], i)
		
		sampleWrite = open(sampleOut, 'w')
		sampleWrite.write("FAMID INDID FAT_ID MAT_ID SEX PHENO\n")
		for i in callData:
			headerCheck = i
			if headerCheck[0:11] == 'probeset_id':
				headerArray = headerCheck.split()
				for j in range(1, len(headerArray)):
					#idArray = headerArray[j].split('_')
					#id = idArray[4].replace(".CEL", '')
					id = headerArray[j]
					if id in sexDic:
						#sampleWrite.write(id + ' ' + id + " 0 0 " + sexDic[id] + " 1\n")
						sampleWrite.write(sexDic[id] + "\n")
		sampleWrite.close()
		
	
	def annoExtract(self, annoIn, annoOut):	
		print "CLASS: Running annoExtract()..."
		
		annoData = fileIO.fileRead(annoIn)
		annoWrite = open(annoOut, 'w')
		for i in annoData:
			annoText = i
			if annoText[0] != '#':
				annoArray = annoText.replace('"', '').split(',')
				data = (annoArray[0] + '\t' + annoArray[1] + '\t' + annoArray[2] + '\t' + 
					annoArray[4] + '\t' + annoArray[5] + '\t' + annoArray[7] + '\t' + 
					annoArray[9] + '\t' + annoArray[11] + '\t' + annoArray[12] + '\t' + 
					annoArray[13] + '\t' + annoArray[14])
				annoWrite.write(data + '\n')
		annoWrite.close()
		