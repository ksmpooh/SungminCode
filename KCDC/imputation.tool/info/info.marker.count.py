#awk '($2>=0.8){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
#awk '($2>=0.8) && $3>=0.05){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
#awk '($2>=0.8) && ($3<0.05 && $3>=0.01){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
#awk '($2>=0.8) && $3<0.01){print $0}' IMPUTE4.chr1.ID_info_maf.buffer.txt | wc -l
import os,glob

def main():
    dfs = glob.glob("*txt")
    print(dfs)
    for df in dfs:
        print(df)
        os.system("awk \'($2>=0.8){print $0}\' %s | wc -l"%df)
        os.system("awk \'($2>=0.8) && ($3>=0.05){print $0}\' %s | wc -l"%df)
        os.system("awk \'($2>=0.8) && ($3<0.05 && $3>=0.01){print $0}\' %s | wc -l"%df)
        os.system("awk \'($2>=0.8) && ($3<0.01){print $0}\' %s | wc -l"%df)
        

main()
        