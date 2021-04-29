import os,glob
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]

def main():
    outDir = "merge/"
    os.system("mkdir "+outDir)
    for pheno in phenoList:
        #ALT_logz.onlyID.txt
        os.system("cat *%s*.txt > %s%s.ID.txt"%(pheno,outDir,pheno))


main()