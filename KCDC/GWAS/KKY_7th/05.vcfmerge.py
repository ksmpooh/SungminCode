import glob,os

wdir = "/DATA/smkim/KKY/05.imputation/"
#inDir = wdir + "INPUTs/"
inDir = "/BDATA/smkim/JG/05.imputation/INPUTs/"
shDir = wdir + "SCRIPTs/"
shDir = shDir + "05.vcf.merge/"
os.system("mkdir "+shDir)
outDir = wdir + "OUTPUTs/"
vcfDir = outDir + "05.vcf.filter/"
vcfMergeDir = outDir + "05.vcf.merge/"
os.system("mkdir " +vcfMergeDir)
tool = "bcftools"
out_prefix = "KKY.7th.20211125."
def vcf_merge_command(vcflist,chrom):
    out = open(shDir+chrom+".merge.sh","w")
    out.write("bcftools concat")
    outname = vcfMergeDir+out_prefix+chrom+".vcf.gz"
    for vcf in vcflist:
        out.write(" %s"%vcf)
    out.write(" -Oz >%s\n"%(outname))
    out.write("tabix -f -p vcf %s"%outname)
    out.close()


def main():
    print("main..")
    ref_file = open(inDir + "imputation.POS.auto_final_20200303.txt","r")
    ref_list = [x.replace('\n','').split('\t') for x in ref_file]
    ref_file.close()
    #Ref_Pos	Imputation_Pos
    #chr1_10178_5000000	chr1_10177_5000000
#    chroms = ["chr"+str(x) for x in range(1,22+1)]
    count = 0
    chr1 = "chr1"
    vcfmerge = []
    for ref_pos,imp_pos in ref_list[1:]:
        chr2,front,tail = imp_pos.split('_')
        ref = "%s.%s_%s"%(chr2,front,tail)
        if ref == "chr22.51000001_51244237":
            vcfmerge.append(''.join(glob.glob(vcfDir + "*%s*gz"%ref)))
            vcf_merge_command(vcfmerge,chr1)
            print("%s size: %s"%(chr1,str(len(vcfmerge))))
            count = count + len(vcfmerge)
            break
        elif chr1 != chr2:
            vcf_merge_command(vcfmerge,chr1)
            print("%s size: %s"%(chr1,str(len(vcfmerge))))
            count = count + len(vcfmerge)
            chr1 = chr2
            vcfmerge=[]
        ref = "%s.%s_%s"%(chr2,front,tail)
        vcfmerge.append(''.join(glob.glob(vcfDir + "*%s*gz"%ref)))
    print("count : %s"%str(count))


main()

#print("count : %s"%str(count))


