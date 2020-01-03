'''
Created on 2017. 11. 6.

@author: HMY
'''

import os, sys
from fileInOut import fileInOut
from bimClass import bimClass
fileIO = fileInOut()
bimCall = bimClass()


#def convertToBim(inDir, resultDir, annoIn, dataIn, dataOut, plinkTool):
def convertToBim(annoIn, dataIn, dataOut, plinkTool):
	print "Running convertToBim()..."
	
	#annoIn = inDir + "Axiom_KORV1_0.na34.annot.extract.txt"
	#annoIn = inDir + "Axiom_KORV1_1.na35.annot.extract.txt"
	
	print "Ori. bim copy..."
	bimIn = dataIn + ".bim"
	oriIn = dataIn + ".ori.bim"
	os.system("cp " + bimIn + ' ' + oriIn)
	
	print "RAW convert bim..."
	bimOut = dataIn + ".bim"
	bimCall.convertToBim(annoIn, oriIn, bimOut)
	
	os.system(plinkTool + " --bfile " + dataIn + " --make-bed --out " + dataOut)
	os.system("cp " + oriIn + ' ' + bimIn)
		
	
def main():
	print "MAIN: START!!!"
	

	annoIn = sys.argv[1]
	dataIn = sys.argv[2]
	dataOut = sys.argv[3]

	#inDir = "/RDATA9/myhwang/KCHIP_NC/INPUTs/"
	#resultDir = "/DATA/yjkim/COV/PHARM/"
	
	plinkTool = "plink"
	
	#convertToBim(inDir, resultDir, annoIn, dataIn, dataOut, plinkTool)
	convertToBim(annoIn, dataIn, dataOut, plinkTool)
	
	
main()
