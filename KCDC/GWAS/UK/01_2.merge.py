import os,glob
#phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
phenoList = ["GLU_inv", "HbA1c_inv", "DBP.nRES", "SBP.nRES","PP.nRES","MAP.nRES"]
def main():
    outDir = "merge/"
    os.system("mkdir "+outDir)
    for pheno in phenoList:
        #ALT_logz.onlyID.txt
        #no.covariates.header.txt
        #os.system("cat *%s*.txt > %s%s.ID.txt"%(pheno,outDir,pheno))
        os.system("cat *%s*.txt | grep -v \"rsid\" | awk \'{print $1,$3,$4,$21}' > ./merge/UKB_%s_all.MERGE.for.manhattan.txt"%(pheno,pheno))

main()