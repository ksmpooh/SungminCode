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
import os,glob
phenoList = ["GLU_inv", "HbA1c_inv", "DBP.nRES", "SBP.nRES","PP.nRES","MAP.nRES"]
def forMan():
    outDir = "./forManhattanPlot/"
    for pheno in phenoList:
        input =  "UKB_%s_ALL_MERGE.txt"%pheno
        output = outDir + "UKB_%s_ALL_MERGE_for.manhattan.txt"%pheno
        os.system("grep -v \"rsid\" %s | awk \'$9 >= 0.8 && $19 >= 0.01{print $1,$3,$4,$21}\' > %s"%(input,output))


forMan()



import os,glob
phenoList = ["GLU_inv", "HbA1c_inv", "DBP.nRES", "SBP.nRES","PP.nRES","MAP.nRES","HT"]
def newMERGE():
    #outDir = "./forManhattanPlot/"
    outDir = "./MERGE/"
    os.system("mkdir %s"%(outDir))
    for pheno in phenoList:
        print("pheno : %s"%pheno)
        header = "header.txt"
        output = outDir + "UKB_%s_ALL_MERGE.txt"%pheno
        os.system("cp %s %s"%(header,output))
        os.system("cat ./%s/*%s*MERGE.txt | grep -v \"rsid\" >> %s"%(pheno,pheno,output))
        output = outDir + "UKB_%s_ALL_MERGE_filINFO0.8.txt"%pheno
        os.system("cp %s %s"%(header,output))
        os.system("cat ./%s/*%s*MERGE_filINFO0.8.txt | grep -v \"rsid\" >> %s"%(pheno,pheno,output))


newMERGE()

#  1            2       3         4         5       6   7       8                              9    10          11          12          13          14      15      16      17          18      19      20                      21                      22                  23                      24                  25      26          27
#alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total all_maf missing_data_proportion frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1 comment	bgenINFO	bgenMAF
import os,glob
phenoList = ["GLU_inv", "HbA1c_inv", "DBP.nRES", "SBP.nRES","PP.nRES","MAP.nRES","HT"]
def forMan2():
    outDir = "./forManhattanPlot/"
    os.system("mkdir %s"%outDir)
    for pheno in phenoList:
        input =  "UKB_%s_ALL_MERGE_filINFO0.8.txt"%pheno
        output = "UKB_%s_ALL_MERGE_filINFO0.8_MAF0.01.txt"%pheno
        header = "../header.txt"
        os.system("cp %s %s"%(header,output))
        os.system("awk '$27>=0.01{print $0}' %s >> %s"%(input,output))
        output2 = outDir + "UKB_%s_ALL_MERGE_filINFO0.8_MAF0.01.for.manhattan.txt"%pheno
        #os.system("grep -v \"rsid\" %s | awk \'$9 >= 0.8 && $19 >= 0.01{print $1,$3,$4,$21}\' > %s"%(input,output))
        os.system("grep -v \"rsid\" %s | awk \'{print $1,$3,$4,$21}\' > %s"%(output,output2))


forMan2()