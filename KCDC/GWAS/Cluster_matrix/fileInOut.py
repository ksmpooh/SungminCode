'''
Created on 2013. 10. 18.

@author: aldud20011



Code study..smkim on 2019.07.10
'''

import os, gzip

class fileInOut():
	#print "CLASS: fileInOut"
	
	def makeDir(self, outDir):   
		#print "METHOD: makeDir()..."
		
		os.system("mkdir " + outDir)
	
	def fileRead(self, fileIn):
		#print "METHOD: fileRead()..."
		
		fileList = open(fileIn, 'r')
		inData = [r.replace('\r', '').replace('\n', '') for r in fileList]
		
		return inData
	
	def gzRead(self, fileIn):
		print "METHOD: gzRead()..."
		
		fileList = gzip.open(fileIn, 'r')
		inData = [r.replace('\r', '').replace('\n', '') for r in fileList]
		
		return inData
		
	def fileWrite(self, writeData, fileOut):
		print "METHOD: fileWrite()..."
		
		outFile = open(fileOut, 'w')
		
		for i in writeData:
			outFile.write(i + '\n')
		outFile.close()
	
	
	def makeTemp(self, tempDir, dataType):
		print "METHOD: makeTemp()..."
		
		tempDir = tempDir + dataType + '/'
		self.makeDir(tempDir)
		
		return tempDir
	