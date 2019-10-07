'''
Created on 2018. 7. 25.

@author: myhwang
'''

import os, sys
from fileInOut import fileInOut
fileIO = fileInOut()

#def indelChange(inDir, sampleDir, plinkTool):
def indelChange(indelIn, bimIn, bimOut, plinkTool):
	print "Running indelChange..."
	
	#indelIn = inDir + "Axiom_KORV1_0.na34.annot.extract.onlyINDEL.txt"
	#bimIn = sampleDir + "KNIH.RAW.35K.2nd.QCed.changeID"
	#bimOut = sampleDir + "KNIH.RAW.35K.2nd.QCed.changeID.indel"
	
	os.system("cp " + bimIn + ".bim " + bimIn + ".ori.bim")
	
	indelDic = {}
	indelData = fileIO.fileRead(indelIn)
	for i in indelData[1:]:
		indelArray = i.split()
		indelDic.setdefault(indelArray[0], '\t'.join(indelArray[1:8]))
	
	bimData = fileIO.fileRead(bimIn + ".bim")
	bimWrite = open(bimIn + ".bim", 'w')
	for i in bimData:
		bimArray = i.split()
		bimCHR = bimArray[0]
		bimID = bimArray[1]
		bimPOS = bimArray[3]
		bimA = bimArray[4]
		bimB = bimArray[5]
		
		if bimID in indelDic:
			indelArray = indelDic[bimID].split()
			oriB = indelArray[1]
			indelPOS = indelArray[2]
			indelA = indelArray[3]
			indelB = indelArray[4]
			if bimA == '-' and bimB == '0':
				bimA = indelA
				bimPOS = indelPOS
			elif bimA == '0' and bimB == '-':
				bimB = indelA
				bimPOS = indelPOS
			elif bimA == '0' and bimB == oriB:
				bimB = indelB
				bimPOS = indelPOS
			elif bimA == oriB and bimB == '0':
				bimA = indelB
				bimPOS = indelPOS	
			elif bimA == '-':
				bimA = indelA
				bimB = indelB
				bimPOS = indelPOS
			elif bimB == '-':
				bimB = indelA
				bimA = indelB
				bimPOS = indelPOS				
				
		bimWrite.write(bimCHR + '\t' + bimID + "\t0\t" + bimPOS + '\t' + bimA + '\t' + bimB + '\n')
	bimWrite.close()

	os.system(plinkTool + " --bfile " + bimIn + " --make-bed --out " + bimOut)
	os.system("cp " + bimIn + ".ori.bim " + bimIn + ".bim")



def main():
	print "START: main"
	
	#inDir = "/jdata/scratch/myhwang/KCHIP_35K_2nd/INPUTs/"
	#sampleDir = "/jdata/scratch/myhwang/KCHIP_35K_2nd/QCed_35K/"
	
	indelIn = sys.argv[1]
	bimIn = sys.argv[2]
	bimOut = sys.argv[3]
	
	
	plinkTool = "plink"
	#indelChange(inDir, sampleDir, plinkTool)
	indelChange(indelIn, bimIn, bimOut, plinkTool)
	
main()






