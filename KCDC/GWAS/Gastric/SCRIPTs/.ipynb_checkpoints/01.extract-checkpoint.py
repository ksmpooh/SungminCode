##extract each chr

import os
chr_input_dir = "/DATA/KCHIP_2019/Gastric/OAS_Plink/"
def extract():
	for i in range(1,22+1):
		os.system("plink --bfile 2nd_New_merge/NEW_MERGE_rmKing --extract "+chr_input_dir+"chr"+str(i)+".txt --make-bed --out 2nd_New_merge/NEW_MERGE_rmKing.chr"+str(i))
extract()




##1000genome DATA loc : /DATA/KCHIP_2019/Gastric/OAS_Plink

### to make intersect ID 1000genome with new_merge


gm_input_dir="/DATA/KCHIP_2019/Gastric/OAS_Plink/"
merge_input_dir="/DATA/smkim/Gastric/QC_2nd/OUTPUTs/2nd_New_merge/"
out_dir = "/DATA/smkim/Gastric/QC_2nd/INPUTs/Merge_with1000genomeID/"
def extract_intersect_ID():
	for i in range(1,22+1):
		os.system("awk '{print $2}' "+gm_input_dir +"ALL.chr"+str(i)+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bim > "+merge_input_dir+"1000genomeSNP_ID_chr"+str(i)+".txt")
		os.system("awk '{print $2}' "+merge_input_dir +"NEW_MERGE_rmKing.chr"+str(i)+".bim > "+merge_input_dir+"mergeSNP_ID_chr"+str(i)+".txt")
		os.system("cat "+merge_input_dir+"1000genomeSNP_ID_chr"+str(i)+".txt "+merge_input_dir+"mergeSNP_ID_chr"+str(i)+".txt |sort|uniq -c|awk '$1 == 2{print $2}'  > "+merge_input_dir+"/chr"+str(i)+"_intersectSNP_ID.txt") 
extract_intersect_ID()


def plink_extract_intersect_ID():
	for i in range(1,22+1):
		os.system("plink --bfile "+merge_input_dir+"NEW_MERGE_rmKing.chr"+str(i)+" --extract "+merge_input_dir+"chr"+str(i)+"_intersectSNP_ID.txt --make-bed --out "+merge_input_dir+"NEW_MERGE_rmKing.chr"+str(i)+".intersect")
		os.system("plink --bfile "+gm_input_dir +"ALL.chr"+str(i)+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --extract "+merge_input_dir+"chr"+str(i)+"_intersectSNP_ID.txt --make-bed --out "+merge_input_dir+"All.chr"+str(i)+".phase3.intersect")
plink_extract_intersect_ID()



