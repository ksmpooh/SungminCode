import os
from fileInOut import fileInOut
fileIO = fileInOut()

inDir = "/DATA/smkim/Gastric/last_qc_2nd/INPUTs/"
outDir = "/DATA/smkim/Gastric/last_qc_2nd/OUTPUTs/"

case = "CASE/"
control = "CONTROL/"






def extract_kingID():
	
	missData = fileIO.fileRead(outDir+case+"KNIH.RAW.Gastric.2nd_rmSNP_rmSample_rmaffy_rmdup_convert_indel_flip_miss.imiss")
	missArray = [j.split() for j in missData[1:]]
	missDic = {j[1]:j[3] for j in missArray}

	with open(outDir+case+"case_kingID_selection.txt", 'w') as kingWrite:
		kingList = fileIO.fileRead(outDir+case+"KNIH.RAW.Gastric.2nd_rmSNP_rmSample_rmaffy_rmdup_convert_indel_flip.kin0")
		for j in kingList[1:]:
			kingArray = j.split()
                        kingID1 = kingArray[1]
                        kingID2 = kingArray[3]

                        kingM1 = int(missDic[kingID1])
                        kingM2 = int(missDic[kingID2])

                        if kingM1 < kingM2:
                        	kingWrite.write(kingID2+"\t"+kingID2 + '\n')
                        else:
                                kingWrite.write(kingID1 +"\t"+kingID1 +  '\n')
extract_kingID()



