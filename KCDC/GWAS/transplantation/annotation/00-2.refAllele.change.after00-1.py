# 20220316

import os,glob,sys

wDir = "/BDATA/smkim/JG/anno/OUTPUTs/merge/"
inDir = wDir + "ori_withRef/"
outDir = wDir + "afterCheckRefAllele/"

#16	60291	16:60291_T/C	T	C	.	.	.	T
#16	60375	16:60375_A/AC	A	AC	.	.	.	A
#16	60842	16:60842_A/C	A	C	.	.	.	A
#16	61349	16:61349_G/A	G	A	.	.	.	G
#16	61730	16:61730_C/G	C	G	.	.	.	C
#16	61977	16:61977_C/G	C	G	.	.	.	C


##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
def main():
    vcfs = glob.glob(inDir + "*.txt")
        
    for vcf in vcfs:
        print(vcf)
        outPath = vcf.replace(inDir,outDir).replace("_withREF.txt","_checkREFallele.vcf")
        inData = open(vcf,"r")
        outData = open(outPath,"w")
        #outData.write(inData.readline())
        #outData.write(inData.readline())
        outData.write("##fileformat=VCFv4.1\n")
        outData.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        while 1:
            line = inData.readline()
            if not line: 
                break
            line = line.replace("\n","")
            chrom,pos,ID,a1,a2,a,b,c,ref = line.split("\t")
            if a1 == ref:
                outData.write("%s\t%s\t%s\t%s\t%s\t.\t.\t.\n"%(chrom,pos,ID,a1,a2))
            else:
                outData.write("%s\t%s\t%s\t%s\t%s\t.\t.\t.\n"%(chrom,pos,ID,a2,a1))

main()