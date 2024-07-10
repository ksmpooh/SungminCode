#######AD MC extract
# python3 script.py [input TRGT VCF]

import os,glob,sys
import re
from pathlib import Path


def mk_df_out(out,mc,ap,am):
    #df = open(,"r")
    mc_in = open(mc,"r")
    ap_in = open(ap,"r")
    am_in = open(am,"r")
    out = open(out,"w")
    
    count = 0
    while 1:
        count = count + 1
        mc_line = mc_in.readline().strip()
        ap_line = ap_in.readline().strip()
        am_line = am_in.readline().strip()
        if not mc_line:
            break
        mc_tmp = mc_line.split()
        ap_tmp = ap_line.split()
        am_tmp = am_line.split()
        if count == 1:
            header = 'ID\tCHROM\tPOS\tTRID\tMOTIFS\tSTRUC\tALT\tAllele\tMC\tAP\tAM\n' # alt "." 경우 STR가 없는 거
            sampleID = mc_tmp[-1].split("=")[0]
            out.write(header)
        mc_allele = mc_tmp[-1].split("=")[1].split(",")
        ap_allele = ap_tmp[-1].split("=")[1].split(",")
        am_allele = am_tmp[-1].split("=")[1].split(",")
        if mc_allele.count(".") == 0:
            #print("%s\t%s\t%s\t%s\t%s\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_1",mc_allele[0],ap_allele[0],am_allele[0]))
            out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_1",mc_allele[0],ap_allele[0],am_allele[0]))
            out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_2",mc_allele[1],ap_allele[1],am_allele[1]))
        else:
            out.write("%s\t%s\t%s\tNA\tNA\tNA\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_1"))
            out.write("%s\t%s\t%s\tNA\tNA\tNA\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_2"))
        '''
        if mc_allele.count(".") == 0:
            out.write("%s\t%s\t%s\t%s\t%s\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_1",mc_allele[0],ap_allele[0]))
            out.write("%s\t%s\t%s\t%s\t%s\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_2",mc_allele[1],ap_allele[1]))
        else:
            out.write("%s\t%s\t%s\tNA\tNA\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_1"))
            out.write("%s\t%s\t%s\tNA\tNA\n"%(sampleID,"\t".join(mc_tmp[:-1]),"allele_2"))
        '''
    out.close()

#chr1:31555:chr1_31555_31570:AAAAT:(AAAAT)n
#chr1    35488   chr1_35488_35504        chr1_35488_35504        AAAT    <STR3>,<STR1>   NIH20N2000078=0/0       NIH20N2038392=0/0       NIH20N2042469=0/0       NIH20N2052743=0/0


def main():
    #inVCF_long = wDir + "Revio.STR.pbmm2_hg38_withunmapped_trgt_genotype_gangstr.sorted.merged_withall.vcf.gz"
    #inVCF_short = wDir + "Shortread.STR.goast_hg38_withunmapped_EH_genotype_gangstr.merged_withall.vcf.gz"
    inVCF = sys.argv[1]
    outDir = "./Quality_check/"
    if Path(outDir).exists() == False:
        os.system("mkdir %s"%outDir)
    


    short_out = outDir + inVCF.replace(".vcf.gz",".MC_AP_AM.txt")
    #short_tmp = inVCF_short.replace(".vcf.gz",".tmp")

    short_mc = outDir + inVCF.replace(".vcf.gz",".tmp.mc")
    short_ap = outDir + inVCF.replace(".vcf.gz",".tmp.ap")
    short_am = outDir + inVCF.replace(".vcf.gz",".tmp.am")

    os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%VARID\\t%%REPID\\t%%RU\\t%%ALT\\t%%FILTER[\\t%%SAMPLE=%%GT]\\n' %s > %s"%(inVCF,short_mc))
    os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%TRID\\t%%MOTIFS\\t%%STRUC\\t%%ALT[\\t%%SAMPLE=%%MC]\\n' %s > %s"%(inVCF,short_mc))
    os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%TRID\\t%%MOTIFS\\t%%STRUC\\t%%ALT[\\t%%SAMPLE=%%AP]\\n' %s > %s"%(inVCF,short_ap))
    os.system("bcftools query -f '%%CHROM\\t%%POS\\t%%TRID\\t%%MOTIFS\\t%%STRUC\\t%%ALT[\\t%%SAMPLE=%%AM]\\n' %s > %s"%(inVCF,short_am))
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%SO\t%ADSP\t%ADFL\t%ADIR]\n'

    ## make df
    print("make df......")
    mk_df_out(long_out,long_mc,long_ap,long_am)



main()
