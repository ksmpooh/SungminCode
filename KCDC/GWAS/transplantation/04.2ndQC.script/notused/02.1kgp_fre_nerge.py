import os
from fileInOut import fileInOut
fileIO = fileInOut()

POP_list =  ['afrID.samples', 'amrID.samples',  'easID.samples' , 'eurID.samples', 'sasID.samples']
inDir = "/DATA/smkim/Gastric/last_qc_2nd/INPUTs/"
outDir = "/DATA/smkim/Gastric/last_qc_2nd/OUTPUTs/"

case = "CASE/"
control = "CONTROL/"
kgp = "1KGP/"


OASDir = "/DATA/KCHIP_2019/Gastric/OAS_Plink/"

prefix = ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes"

#plink duplicate  for each chr
def rm_duplicate_each_pop() :
	for i in range(1,22+1):
		i = str(i)
		os.system("plink --bfile "+OASDir+"ALL.chr"+i+prefix+" --missing --freq --a1-allele "+inDir+kgp+"ALL"+prefix+".REF.txt --list-duplicate-vars --out "+inDir+kgp+"tmp")

		missData = fileIO.fileRead(inDir+kgp+"tmp.lmiss")
		missArray = [j.split() for j in missData[1:]]
		missDic = {j[1]:j[4] for j in missArray}

		with open(inDir+kgp+"chr"+i+"_duplicateID.txt", 'w') as dupWrite:
			dupList = fileIO.fileRead(inDir+kgp+"tmp.dupvar")
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


#		os.system("awk '{print $4}' "+inDir+kgp+"tmp.dupvar > "+inDir+kgp+"chr"+i+"_duplicateID.txt")
		os.system("awk '$5 < 0.05{print $2}' "+inDir+kgp+"tmp.frq > "+inDir+kgp+"chr"+i+"_frqID.txt")
		os.system("cat "+inDir+kgp+"chr"+i+"_duplicateID.txt "+inDir+kgp+"chr"+i+"_frqID.txt > "+inDir+kgp+"1kgp_rmMARKER_chr"+i+".txt")
		os.system("plink --bfile "+OASDir+"ALL.chr"+i+prefix+" --maf 0.05 --exclude "+inDir+kgp+"1kgp_rmMARKER_chr"+i+".txt --make-bed --out "+outDir+kgp+"1KGP_chr"+i) 

#os.system("mkdir "+inDir + "pop_merge_list/")
#os.system("mkdir /DATA/smkim/Gastric/reQC_2nd/SCRIPTs/pop_merge/")

def merge():
	out = open(inDir+kgp + "1kgp_merge_list.txt","w")
	for i in range(2,22+1):
		i = str(i)
		out.write(outDir+kgp+"1KGP_chr"+i+".bed\t"+outDir+kgp+"1KGP_chr"+i+".bim\t"+outDir+kgp+"1KGP_chr"+i+".fam\n")
	out.close()
	os.system("plink --bfile "+outDir+kgp+"1KGP_chr1 --merge-list "+inDir+kgp+"1kgp_merge_list.txt --allow-no-sex --make-bed --out "+outDir+kgp+"1kgp_merge")


rm_duplicate_each_pop()
merge()


