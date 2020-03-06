import os,sys
from fileInOut import fileInOut
fileIO = fileInOut()







def extract_duplicate_ID(file):
#	os.system("plink --bfile "+outDir+case+"KNIH.RAW.Gastric.2nd_rmSNP_rmSample_rmaffy --missing --list-duplicate-vars --out "+outDir+case+"KNIH.RAW.Gastric.2nd_rmSNP_rmSample_rmaffy")
        os.system("plink --bfile "+file+" --missing --list-duplicate-vars --out "+file+".miss-dup")

	missData = fileIO.fileRead(file+".miss-dup.lmiss")
	missArray = [j.split() for j in missData[1:]]
	missDic = {j[1]:j[4] for j in missArray}

#	with open(outDir+case+"case_duplicateSNPID.txt", 'w') as dupWrite:
        with open("case_duplicateSNPID.txt", 'w') as dupWrite:

		dupList = fileIO.fileRead(file+".miss-dup.dupvar")
		for j in dupList[1:]:
			dupArray = j.split()
                        dupID1 = dupArray[3]
                        dupID2 = dupArray[4]

                        dupM1 = float(missDic[dupID1])
                        dupM2 = float(missDic[dupID2])

                        if dupM1 < dupM2:
                        	dupWrite.write(dupID2 + '\n')
                        else:
                                dupWrite.write(dupID1 + '\n')
def main():
	name = sys.argv[1]
	extract_duplicate_ID(name)

main()
