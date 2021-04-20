#awk '($2>=0.8){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
#awk '($2>=0.8) && $3>=0.05){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
#awk '($2>=0.8) && ($3<0.05 && $3>=0.01){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
#awk '($2>=0.8) && $3<0.01){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
#impute5.chr1.ID_info_maf.eagle_vcftovcf_pasing.usingKRG1KGP.withMAF.txt', 
#'impute5.chr1.ID_info_maf.eagle_plinktohap_phasing.usingKRG1KGP.withMAF.txt', 
# 'impute5.chr1.ID_info_maf.shapeit4_pasing.usingKRG1KGP.withMAF.txt', 
# 'minimac4.eagle_phasing.default_window.KRG1KGP.withMAF.txt', 
# 'minimac4.eagle_phasing.use_window.KRG1KGP.withMAF.txt', 
# 'minimac4.shapeit4_phasing.use_window.KRG1KGP.withMAF.txt', 
# 'minimac4.shapeit4_phasing.Defualt_window.KRG1KGP.withMAF.txt', 
# 'IMPUTE4.chr1.ID_info_maf.buffer.withMAF.tx
import os,glob

def main():
    dfs = glob.glob("*txt")
    print(dfs)
    for df in dfs:
        print(df)
        os.system("cat %s |grep -v ID | wc -l"%df)
        os.system("awk \'($3>=0.05){print $0}\' %s | grep -v ID|wc -l"%df)
        os.system("awk \'($3<0.05 && $3>=0.01){print $0}\' %s | wc -l"%df)
        os.system("awk \'($3<0.01){print $0}\' %s | grep -v ID|wc -l"%df)
        os.system("awk \'($2>=0.8){print $0}\' %s | grep -v ID|wc -l"%df)
        os.system("awk \'($2>=0.8) && ($3>=0.05){print $0}\' %s |grep -v ID| wc -l"%df)
        os.system("awk \'($2>=0.8) && ($3<0.05 && $3>=0.01){print $0}\' %s | grep -v ID | wc -l"%df)
        os.system("awk \'($2>=0.8) && ($3<0.01){print $0}\' %s | grep -v ID | wc -l"%df)
        

main()
        