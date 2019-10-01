import os
POP_list =  ['afrID.samples', 'amrID.samples',  'easID.samples' , 'eurID.samples', 'sasID.samples']
inDir = "/DATA/smkim/Gastric/QC_2nd/INPUTs/"
outDir = "/DATA/smkim/Gastric/QC_2nd/OUTPUTs/Last_merge/"

for pop in POP_list:
        os.system("plink --bfile "+outDir+"1Kgenome_"+pop+" --maf 0.05 --freq --a1-allele "+inDir+"change_SNPID_type_Axiom_KOR.annot.extract.addINDEL.FINAL.REF_onlyOne.txt --make-bed --out "+outDir+"1Kgenome_"+pop+"_frq")

