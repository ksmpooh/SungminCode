import os
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
def main():
    for pheno in phenoList:
        #UKB_130K_GGT_logz_ALL.merge_20210429.txt
        input = "UKB_130K_%s_ALL.merge_20210429.txt"%pheno
        output = "UKB_130K_%s_marker.ID_130K.txt"%pheno
        #grep -v alternate_ids
        os.system("grep -v alternate_ids %s | cut -d\" \" -f 1 > %s"%(input,output))

main()        

    
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

import os
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
def main():
    
    for pheno in phenoList:
        input1 = "/BDATA/myhwang/KBA_130K/11_UKB/RESULTs/assoMERGE_130K_all/chr_pos/UKB_130K_%s_ALL.merge_20210429_chr.pos.ID.txt"%pheno
        input2 = "/BDATA/smkim/UK/130K/marker_ori/onlyID/chr_pos/%s.ID.txt"%pheno
        output = "/BDATA/smkim/UK/130K/marker.check/%s_markercheck.txt"%pheno
        os.system("cat %s %s |grep -v \"chr\" | sort | uniq -c| awk '$1==1{print $2}' > %s"%(input1,input2,output))

main()


####################
#alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total all_maf missing_data_proportion frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1 comment
#$3 $4

import os,glob
phenoList = ["ALT_logz", "GLU_inv", "HbA1c_inv", "AST_logz", "GGT_logz"]
def main():
    for pheno in phenoList:
        input = "UKB_130K_%s_ALL.merge_20210429.txt"%pheno
        #input = glob.glob("*%s*"%pheno)
        output = "chr_pos/UKB_130K_%s_ALL.merge_20210429_chr.pos.ID.txt"%pheno
        #os.system("awk '{print $3\":\"$4}' %s > %s"%(input,output))
        os.system("cut -d\" \" -f3-4 %s | awk '{print $1\":\"$2}' > %s"%(input,output))

main()