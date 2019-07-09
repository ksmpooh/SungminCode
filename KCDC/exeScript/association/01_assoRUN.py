'''
Created on 2019. 5. 31.

@author: myhwang
'''

import glob
from fileInOut import fileInOut
from assoClass import assoClass
fileIO = fileInOut()
assoCall = assoClass()


def	qtRUN(inDir, scriptDir, filDir, assoDir, epactsTool):
	print "Running qtRUN..."
	
	phenoList = ["CPD"]
	for i in phenoList:
		phenoType = i
		
		shDir = fileIO.makeTemp(scriptDir, "assoQT_" + phenoType)
		traitDir = fileIO.makeTemp(assoDir, phenoType)
		
		pedIn = inDir + "KCHIP_SMK_CPD_201906.ped"
		covType = " --cov SEX --cov AGE --cov AS --cov DS --cov NC"
		
		assoCall.qtRUN(shDir, filDir, traitDir, pedIn, covType, phenoType, epactsTool, chipType="V1")
		assoCall.qtRUN(shDir, filDir, traitDir, pedIn, covType, phenoType, epactsTool, chipType="V2")
			

def assoMERGE(inDir, assoDir, mergeDir):
	print "Running assoMERGE..."
	
	chunkIn = inDir + "imputation.IMPUTE4.POS.50K_20181114_Final.txt"
	
	phenoList = ["CPD"]
	for i in phenoList:
		phenoType = i
		traitDir = fileIO.makeTemp(assoDir, phenoType)
		assoCall.assoMERGE(traitDir, mergeDir, chunkIn, phenoType)
		

def manQQ(scriptDir, mergeDir, plotTool):
	print "Running assoMERGE..."
	
	shDir = fileIO.makeTemp(scriptDir, "manQQ")
	
	phenoList = ["CPD"]
	for i in phenoList:
		phenoType = i
		epactsData = glob.glob(mergeDir + phenoType + "_*.txt")
		epactsList = [r.replace('\r', '').replace('\n', '') for r in epactsData]
		for j in epactsList:
			shOut = j.replace(mergeDir, shDir).replace(".txt", ".sh")
			pdfOut = j.replace(".txt", '')
			with open(shOut, 'w') as shWrite:
				shWrite.write(plotTool + " -in " + j + " -out " + pdfOut + '\n')

		
def main():
	
	#serverType = "mac"
	serverType = "oas"
	if serverType == "mac":
		defaultDir = "/Users/myhwang/workspace/36_JointResearch/src/12_SHJI_Smoking/"
	else:
		defaultDir = "/jdata/scratch/myhwang/Cowork/SHJI/"
		qtDir = "/jdata/scratch/myhwang/KCHIP_136K/03_QT/RESULTs/"
		epactsTool = "/ldata/tools/EPACTS-3.2.6/bin/epacts"
		plotTool = "/ldata/tools/EPACTS-3.2.6/bin/epacts-plot"
	
	
	inDir = defaultDir + "INPUTs/"
	scriptDir = defaultDir + "SCRIPTs/"
	resultDir = defaultDir + "RESULTs/"
	filDir = qtDir + "filterIMPUTE4/"
	
	assoDir = fileIO.makeTemp(resultDir, "assoRESULTs")
	mergeDir = fileIO.makeTemp(resultDir, "assoMERGEs")
	
	#qtRUN(inDir, scriptDir, filDir, assoDir, epactsTool)
	assoMERGE(inDir, assoDir, mergeDir)
	#manQQ(scriptDir, mergeDir, plotTool)
	
main()
