import os, glob, sys
'''
/BDATA/smkim/TOOLs/hla_mapper/bin/linux/hla-mapper dna \
r1=/BDATA/smkim/HLA_seq/00.rawDATA/1stTest/247_S99_L002_R1_001_paired.fastq.gz \
r2=/BDATA/smkim/HLA_seq/00.rawDATA/1stTest/247_S99_L002_R2_001_paired.fastq.gz \
output=/BDATA/smkim/HLA_seq/Tool_test/hla_mapper/shortread/test/ \
sample=NIH19KT0247 db=/BDATA/smkim/TOOLs/hla_mapper/hla-mapper_db_004.1_HLA/ \
threads=32
'''
wDir = "/BDATA/smkim/HLA_seq/"
refDir = wDir + "REF/"
shortread_Dir = wDir + "shortread/"
shortread_raw = wDir +  "00.rawDATA/shortread/"

def HLA_mapper_Shortread(ID_table):

    Tool = "/BDATA/smkim/TOOLs/hla_mapper/bin/linux/hla-mapper"
    db = "/BDATA/smkim/TOOLs/hla_mapper/hla-mapper_db_004.1_HLA/"

    for i in ID_table:
        i = i.split()
        #print(i)
        R1 = shortread_raw + i[2] + "_paired.fastq.gz"
        R2 = R1.replace("_R1_","_R2_")
        outDir = "/BDATA/smkim/HLA_seq/HLA_infer/HLA_mapper/shortread/%s/"%i[0]
        os.system("mkdir %s"%outDir)
        #out = shortread_mapping + "HLA.Shortread.Seq.%s.align.sorted.bam"%i[0]
        #out = "%sHLA.Shortread.Seq.%s.trimmed.%s_align.sorted.bam"%(shortread_mapping,i[0],theme)
        os.system("%s dna r1=%s r2=%s output=%s sample=%s db=%s threads=32"%(Tool,R1,R2,outDir,i[0],db))


def HLA_scan_Shortread():
    print("HLA scan")
    

def main():
    refs = open(wDir + "HLAseq.ID.table_v2.txt","r")
    refs = [s.replace("\n","") for s in refs]
    header = refs.pop(0)

    HLA_mapper_Shortread(refs)