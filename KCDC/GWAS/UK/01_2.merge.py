import os,glob
#phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
phenoList = ["GLU_inv", "HbA1c_inv", "DBP.nRES", "SBP.nRES","PP.nRES","MAP.nRES"]
def main():
    outDir = "MERGE/"
    os.system("mkdir "+outDir)
    for pheno in phenoList:
        #ALT_logz.onlyID.txt
        #no.covariates.header.txt
        #os.system("cat *%s*.txt > %s%s.ID.txt"%(pheno,outDir,pheno))
        #os.system("cat *%s*.txt | grep -v \"rsid\" | awk \'{print $1,$3,$4,$21}' > ./merge/UKB_%s_all.MERGE.for.manhattan.txt"%(pheno,pheno))
        output = outDir + "UKB_%s_ALL_MERGE.txt"%pheno
        header = "no.covariates.header.txt"
        os.system("cp %s %s"%(header,output))
        os.system("cat *%s*.txt | grep -v \"rsid\" >> %s"%(pheno,output))

main()


#ukb57705_imp_v3_s487283_Liver_20210422.FINAL.sample
#ukb57705_imp_v3_s487283_GLY_20210422.FINAL.sample
#ukb57705_imp_v3_s487283_HT_20210517.sample
1                 3           4                                                               9                                                                                                 19
#alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total all_maf missing_data_proportion frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1 comment
def forMan():
    outDir = "forManhattanPlot/"
    for pheno in phenoList:
        input =  "UKB_%s_ALL_MERGE.txt"%pheno
        output = oudtDir + "UKB_%s_ALL_MERGE_for.manhattan.txt"%pheno
        os.system("grep -v \"rsid\" %s | awk \'$9 >= 0.8 & $19 >= 0.01{print $1,$3,$4,$21}\' > %s"%(input,output))
