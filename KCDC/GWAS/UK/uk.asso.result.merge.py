
import os,glob

inDir = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/assoRUN_130K_re/"
outDir = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/assoMERGE_130K_re/"

phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
header = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/ukb.header.txt"

def main():
    for pheno in phenoList:
        make_time = "20210424"
        os.system("cp %s %sGWAS.cal_%s_ALL.merge_%s.txt"%(header,outDir,pheno,make_time))
        os.system("cat %s%s/*.txt | grep -v \"#\" | grep -v \"alternate_ids\" >> %sGWAS.cal_%s_ALL.merge_%s.txt"%(inDir,pheno,outDir,pheno,make_time))


main()