'''
Created on 2016. 1. 13.

@author: HMY
'''

import os, sys
from fileInOut import fileInOut
from prepClass import prepClass
fileIO = fileInOut()
prepCall = prepClass()


def callSplit(inDir, callDir, plinkDir, posSTART, posUnit):
	print "Running callSplit()..."
	
	callIn = callDir + "AxiomGT1.calls.txt"
	callData = fileIO.fileRead(callIn)
	
	mergeOut = inDir + "mergeList.txt"
	posOut = inDir + "posInfo.txt"
	startData, endData = prepCall.posRange(plinkDir, mergeOut, posOut, posSTART, len(callData), posUnit)
	
	startArray = startData.split(',')
	endArray = endData.split(',')
	
	for i in range(0, len(startArray)):
		start = startArray[i]
		end = endArray[i]
		
		callOut = callDir + "AxiomGT1.calls." + str(start) + '_' + str(end) + ".txt"
		if end == "END":
			end = len(callData)
		
		callWrite = open(callOut, 'w')
		for j in range(int(start), int(end)):
			headerCheck = callData[j]
			if headerCheck[0] != '#':
				callWrite.write(headerCheck + '\n')
		callWrite.close()
			

def dataPrep(inDir, callDir, sexIn):
	print "Running dataPrep()..."
	
	posIn = inDir + "posInfo.txt"
	posData = fileIO.fileRead(posIn)
	
	startArray = posData[0].split(',')
	endArray = posData[1].split(',')
	
	callIn = callDir + "AxiomGT1.calls." + str(startArray[0]) + '_' + str(endArray[0]) + ".txt"
	sampleOut = inDir + "KNIH.RAW.info.txt"
	prepCall.sampleID(callIn, sexIn, sampleOut)
	

def scripts(inDir, scriptDir, stepDir, callDir, plinkDir, annoIn):
	print "Running scripts()..."
	
	posIn = inDir + "posInfo.txt"
	posData = fileIO.fileRead(posIn)
	
	startArray = posData[0].split(',')
	endArray = posData[1].split(',')
	
	for i in range(0, len(startArray)):
		start = startArray[i]
		end = endArray[i]
		step1Write = open(stepDir + "CallToPlink_" + str(start) + '_' + str(end) + ".sh", 'w')
		step1Write.write("python " + scriptDir + "CallToPlink.py " + str(start) + ' ' + str(end) +
						' ' + inDir + ' ' + callDir + ' ' + plinkDir + ' ' + annoIn + '\n')
		step1Write.close()
	
	
def main():
	print "MAIN: START!!!"
	
	
	##Run: python DataPrep.py 0 50000 Sample.Info.txt Axiom_KORV1_1.na35.annot.extract.txt	outDir
	##Make Sample.Info.txt
	##FAMID	INDID	FAT_ID	MAT_ID	SEX	PHENO
	##ID001	ID001	0	0	1	1
	##ID002	ID002	0	0	2	1
	
	
	posSTART = int(sys.argv[1])
	posUnit = int(sys.argv[2])
	sexIn = sys.argv[3]
	annoIn = sys.argv[4]
	callDir = sys.argv[5]
	
	scriptDir = ""
	inDir = callDir + "INPUTs/"
	plinkDir = callDir + "PLINK/"
	stepDir = scriptDir + "runSH/"
	
	fileIO.makeDir(inDir)
	fileIO.makeDir(plinkDir)
	fileIO.makeDir(stepDir)
	
	callSplit(inDir, callDir, plinkDir, posSTART, posUnit)
	dataPrep(inDir, callDir, sexIn)
	scripts(inDir, scriptDir, stepDir, callDir, plinkDir, annoIn)
	
main()
