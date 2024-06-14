#######compare
# python3 script.py [input VCF] [eh/trgt] [db] 

import os,glob,sys
import re


wDir = "/BDATA/smkim/STR/01.compare/00.rawDATA/"
wDir = "/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/"

#chr1	98618	chr1_98618_98630	98630	TGAC	(TGAC)n	NIH23F1013274=3,3	NIH23F1093518=3,3	NIH23F1110753=3,3	NIH23F1135991=3,3	NIH23F1140738=3,3
#def main
#CHROM
def remove_chars(text, chars_to_remove):
    for char in chars_to_remove:
        text = text.replace(char, "")
    #print(text)
    return text
#os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%TRID\\t%%MOTIFS\\t%%STRUC\\t%ALT[\\t%%SAMPLE=%%MC]\\n' %s > %s"%(inVCF_long,long_tmp))
def mk_df_out(intmp):
    df = open(intmp,"r")
    out = open(intmp.replace(".tmp",".mkdf.txt"),"w")
    count = 0
    while 1:
        count = count + 1
        line = df.readline().strip()
        if not line:
            break
        tmp = line.split()
        if count == 1:
            header = 'CHROM\tPOS\tTRID\tMOTIFS\tSTRUC\tALT\t' # alt "." 경우 STR가 없는 거
            header = header + "\t".join([s.split("=")[0] for s in tmp[6:]])
            header = header + "\n"
            out.write(header)
        id = "\t".join(tmp[0:5+1])
        out.write("%s\t"%id)
        #print(tmp)
        #print("tmp len : %s"%len(tmp))
        add_line = "\t".join([s.split("=")[1] for s in tmp[6:]])
        out.write("%s\n"%add_line)
    out.close()

#chr1:31555:chr1_31555_31570:AAAAT:(AAAAT)n
#chr1    35488   chr1_35488_35504        chr1_35488_35504        AAAT    <STR3>,<STR1>   NIH20N2000078=0/0       NIH20N2038392=0/0       NIH20N2042469=0/0       NIH20N2052743=0/0

#os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%VARID\\t%%REPID\\t%%RU\\t%%ALT[\\t%%SAMPLE=%%GT]\\n' %s > %s"%(inVCF_short,short_tmp))
def mk_df_out_forEH(intmp):
    df = open(intmp,"r")
    out = open(intmp.replace(".tmp",".mkdf.txt"),"w")
    count = 0
    while 1:
        count = count + 1
        line = df.readline().strip()
        if not line:
            break
        tmp = line.split()
        if count == 1:
            header = 'CHROM\tPOS\tVARID\tREPID\tRU\tALT\tFILTER\t' # alt "." 경우 STR가 없는 거
            header = header + "\t".join([s.split("=")[0] for s in tmp[7:]])
            header = header + "\n"
            out.write(header)
        id = "\t".join(tmp[0:6+1])
        out.write("%s"%id)
        #print(tmp)
        #print("tmp len : %s"%len(tmp))
        #str_n = ("0,"+ keep_numbers_and_delimiter(tmp[5])).split(",")
        if tmp[5] == ".":
            str_n = ["0"]
        else:    
            str_n = ("0,"+ remove_chars(tmp[5],["<", ">", "STR"])).split(",")
        #print(str_n)
        #out.write("%s\t"%str_n)
        for i in tmp[7:]:
            indexing = i.split("=")[1].split("/")
            #print(indexing)
            if indexing[0] == ".":
                a1 = "."
            else:
                a1 = str_n[int(indexing[0])]
            
            if indexing[1] == ".":
                a2 = "."
            else:
                a2 = str_n[int(indexing[1])]
            out.write("\t%s,%s"%(a1,a2))
            #out.write("\t%s,%s"%(str_n[int(indexing[0])],str_n[int(indexing[1])]))
        out.write("\n")
    out.close()




def main():
    inVCF_long = wDir + "Revio.STR.pbmm2_hg38_withunmapped_trgt_genotype_gangstr.sorted.merged_withall.vcf.gz"
    inVCF_short = wDir + "Shortread.STR.goast_hg38_withunmapped_EH_genotype_gangstr.merged_withall.vcf.gz"
    
    long_tmp = inVCF_long.replace(".vcf.gz",".tmp")
    short_tmp = inVCF_short.replace(".vcf.gz",".tmp")
    
    #os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%TRID\\t%%MOTIFS\\t%%STRUC\\t%%ALT[\\t%%SAMPLE=%%MC]\\n' %s > %s"%(inVCF_long,long_tmp))
    os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%VARID\\t%%REPID\\t%%RU\\t%%ALT\\t%%FILTER[\\t%%SAMPLE=%%GT]\\n' %s > %s"%(inVCF_short,short_tmp))

    ## make df
    print("make df......")
    #mk_df_out(long_tmp)
    mk_df_out_forEH(short_tmp)



main()
