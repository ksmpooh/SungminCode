import os
POP_list =  ['afrID.samples', 'amrID.samples',  'easID.samples' , 'eurID.samples', 'sasID.samples']
inDir = "/DATA/smkim/Gastric/QC_2nd/INPUTs/"


for pop in POP_list:
	os.system("plink --bfile /DATA/smkim/Gastric/QC_2nd/OUTPUTs/Last_merge/1kg_merge --keep "+inDir+pop+".txt --make-bed --out /DATA/smkim/Gastric/QC_2nd/OUTPUTs/Last_merge/1Kgenome_"+pop)





