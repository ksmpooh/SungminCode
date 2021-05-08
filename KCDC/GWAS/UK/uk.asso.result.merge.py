
import os,glob

inDir = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/assoRUN_130K_re/"
outDir = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/assoMERGE_130K_re/"

phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
header = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/ukb.header.txt"

def main():
    for pheno in phenoList:
        make_time = "20210429"
        os.system("cp %s %s130K_%s_ALL.merge_%s.txt"%(header,outDir,pheno,make_time))
        os.system("cat %s%s/*.txt | grep -v \"#\" | grep -v \"alternate_ids\" >> %s130K_%s_ALL.merge_%s.txt"%(inDir,pheno,outDir,pheno,make_time))

#main()

def main():
    Dir = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/"
    Dir_1st = Dir + "assoMERGE_130K/"
    Dir_2nd = Dir + "assoMERGE_130K_re/"
    outDir = Dir + "assoMERGE_130K_all/"
    os.system("mkdir "+outDir)
    for pheno in phenoList:
        make_time = "20210429"
        output = "%sUKB_130K_%s_ALL.merge_%s.txt"%(outDir,pheno,make_time)
        input1 = Dir_1st + "UKB_%s_ALL_MERGE_130K.txt"%pheno 
        #UKB_GGT_logz_ALL_MERGE_130K.txt
        input2 = Dir_2nd + "130K_%s_ALL.merge_20210429.txt"%pheno
        #130K_AST_logz_ALL.merge_20210429.txt
        
        os.system("cp %s %s"%(header,output))
        os.system("cat %s %s |  grep -v \"alternate_ids\" >> %s"%(input1, input2, output))

main()



