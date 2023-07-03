'''
/BDATA/smkim/TOOLs/hla_scan/hla_scan_r_v2.1.4 \
-b /BDATA/smkim/HLA_seq/00.rawDATA/1stTest/hla-mapper/NIH19KT0247.unique.bam \
-d /BDATA/smkim/TOOLs/hla_scan/db/HLA-ALL.IMGT \
-v 38 \
-g HLA-DPA1
'''


# python3 HLA-scan.py [input.bam]
import os,sys

tool = "/BDATA/smkim/TOOLs/hla_scan/hla_scan_r_v2.1.4"
db = "/BDATA/smkim/TOOLs/hla_scan/db/HLA-ALL.IMGT"
v = "38"


def main():
    bam = sys.argv[1]
    genes = ["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]

    for gene in genes:
        os.system("%s -b %s -d %s -v %s -g %s"%(tool,bam,db,v,gene))


main()
