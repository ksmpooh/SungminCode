import os
POP_list =  ['afrID.samples', 'amrID.samples',  'easID.samples' , 'eurID.samples', 'sasID.samples']
inDir = "/DATA/smkim/Gastric/QC_2nd/INPUTs/"
outDir = "/DATA/smkim/Gastric/QC_2nd/OUTPUTs/Last_merge/"

for pop in POP_list:
	os.system("plink --bfile "+outDir+"1Kgenome_"+pop+"_frq --extract "+ outDir+"intersect_SNPID_sevenPOP.txt --make-bed --out "+outDir+"1Kgenome_"+pop+"_frq_intersect")

os.system("plink --bfile "+outDir+"Gastric_control_frq --extract "+ outDir+"intersect_SNPID_sevenPOP.txt --make-bed --out "+outDir+"Gastric_control_frq_intersect")
os.system("plink --bfile "+outDir+"Gastric_case_frq --extract "+ outDir+"intersect_SNPID_sevenPOP.txt --make-bed --out "+outDir+"Gastric_case_frq_intersect")


New_outDir = "/DATA/smkim/Gastric/QC_2nd/OUTPUTs/PCA_merge/"

outfile = open(New_outDir + "merge_list.txt","w")

for pop in POP_list:
	outfile.write(outDir+"1Kgenome_"+pop+"_frq_intersect.bed\t"+outDir+"1Kgenome_"+pop+"_frq_intersect.bim\t"+outDir+"1Kgenome_"+pop+"_frq_intersect.fam\n")
outfile.write(outDir +"Gastric_case_frq_intersect.bed\t"+outDir +"Gastric_case_frq_intersect.bim\t"+outDir +"Gastric_case_frq_intersect.fam\n")

outfile.close()












