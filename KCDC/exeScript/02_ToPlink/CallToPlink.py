'''
Created on 2016. 1. 13.

@author: HMY
'''

import sys
from fileInOut import fileInOut
from convertClass import convertClass
fileIO = fileInOut()
convertCall = convertClass()


def convertToMAP(callDir, plinkDir, annoIn, start, end):
	print "Running convertToMAP()..."
	
	
	callName = "AxiomGT1.calls." + str(start) + '_'
	callIn = callDir + "AxiomGT1.calls." + str(start) + '_' + str(end) + ".txt"
	mapOut = plinkDir + "KNIH.RAW." + str(start) + '_' + str(end) + ".map"
	convertCall.covertToMAP(callName, callIn, annoIn, mapOut)
	

def convertToPED(inDir, callDir, plinkDir, start, end):
	print "Running convertToPED()..."
	
	
	callName = "AxiomGT1.calls." + str(start) + '_'
	callIn = callDir + "AxiomGT1.calls." + str(start) + '_' + str(end) + ".txt"
	sampleIn = inDir + "KNIH.RAW.info.txt"
	pedOut = plinkDir + "KNIH.RAW." + str(start) + '_' + str(end) + ".ped"
	convertCall.convertToPED(callName, callIn, sampleIn, pedOut)
	
	
def main():
	print "MAIN: START!!!"
	
	
	start = sys.argv[1]
	end = sys.argv[2]
	inDir = sys.argv[3]
	callDir = sys.argv[4]
	plinkDir = sys.argv[5]
	annoIn = sys.argv[6]
	
	convertToMAP(callDir, plinkDir, annoIn, start, end)
	convertToPED(inDir, callDir, plinkDir, start, end)
	
	
main()
