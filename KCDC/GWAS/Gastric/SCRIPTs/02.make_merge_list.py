import os

indir = "/DATA/smkim/Gastric/QC_2nd/OUTPUTs/2nd_New_merge/"
outdir = "/DATA/smkim/Gastric/QC_2nd/OUTPUTs/"


for file in ["case_control","1kg"]:
	outfile = open(outdir+file+"_merge_list.txt","w")
	if file == "case_control":
		for i in range(2,22+1):
			outfile.write(indir+"NEW_MERGE_rmKing.chr"+str(i)+".intersect.bed\t"+indir+"NEW_MERGE_rmKing.chr"+str(i)+".intersect.bim\t"+indir+"NEW_MERGE_rmKing.chr"+str(i)+".intersect.fam\n")
		outfile.close()
		
	else:
                for i in range(2,22+1):
                        outfile.write(indir+"All.chr"+str(i)+".phase3.intersect.bed\t"+indir+"All.chr"+str(i)+".phase3.intersect.bim\t"+indir+"All.chr"+str(i)+".phase3.intersect.fam\n")
		outfile.close()
os.system("plink --bfile 2nd_New_merge/NEW_MERGE_rmKing.chr1.intersect --merge-list case_control_merge_list.txt --make-bed --allow-no-sex --out Last_merge/Gastric_merge")
os.system("plink --bfile 2nd_New_merge/All.chr1.phase3.intersect --merge-list 1kg_merge_list.txt --make-bed --allow-no-sex --out Last_merge/1kg_merge")





