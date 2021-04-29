import os
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
#UKB_ALT_logz_ALL_MERGE_130K.txt
def main():
    for pheno in phenoList:
        input = "UKB_%s_ALL_MERGE_130K.txt"%pheno
        output = "UKB_%s_marker.ID_130K.txt"%pheno
        #grep -v alternate_ids
        os.system("grep -v alternate_ids %s | cut -d\" \" -f 1 > %s"%(input,output))

#main()        

    
#asso_result  re_marker  target_marke
#GGT_logz.ID.txt
import os
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
def main():
    for pheno in phenoList:
        input1 = "asso_result/UKB_%s_marker.ID_130K.txt"%pheno
        input2 = "target_marker/%s.onlyID.txt"%pheno
        output = "re_marker/UKB_%s_re.marker.ID_130K.txt"%pheno
        os.system("cat %s %s |  sort | uniq -c| awk '$1==1{print $2}' > %s"%(input1,input2,output))

main()