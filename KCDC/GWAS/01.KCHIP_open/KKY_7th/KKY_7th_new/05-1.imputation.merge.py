
import os,glob

minimac4 = "/BDATA/smkim/TOOLs/minimac4"
wdir = "/BDATA/smkim/KKY_7th/"
phasingDir = wdir + "04.phasing/OUTPUTs/02.chr_phasing/"
#phasingDir = "/ADATA/smkim/JG/04.phasing/OUTPUTs/03.5Ksplit/"

imputationDir = wdir + "05.imputation/"


# chr21_39000001_44000000.m3vcf.gz
# genetic_map_chr22_combined_b37_addCHR.txt
#refDir = inDir + "KGP3KRG/"
refDir = "/BDATA/smkim/GWAS/ref/KRG1KGP/m3vcf/"
mapDir = "/BDATA/smkim/GWAS/ref/map/m3vcf/"


vcfDir = imputationDir + "OUTPUTs/01.imputation/"

vcfMergeDir = imputationDir + "OUTPUTs/02.merge/"

os.system("mkdir %s"%vcfMergeDir)

shDir = imputationDir + "SCRIPTs/02.merge/"
os.system("mkdir %s"%shDir)



out_prefix = "KKY.7th.imputation_MINIMAC4."

def vcf_merge_command(vcflist,chrom):
    out = open(shDir+chrom+".merge.sh","w")
    out.write("bcftools concat")
    outname = vcfMergeDir+out_prefix+chrom+".vcf.gz"
    for vcf in vcflist:
        out.write(" %s"%vcf)
    out.write(" -Oz >%s\n"%(outname))
    out.write("tabix -f -p vcf %s"%outname)
    out.close()


#KKY.6th.imputation_MINIMAC4.chr22.16050076_21000000.dose.vcf.gz
def main():
    refs = open(refDir + "ReferencePanel.KRG1KGP.m3vcf.imputation.txt","r")
    refs = [x.replace("\n","") for x in refs]
    chr1 = "chr1"
    count = 0
    vcfmerge = []
    for ref in refs:
        chr2,start,end = ref.split("\t")
        ref = "%s.%s_%s"%(chr,start,end)
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
        ref = "%s.%s_%s"%(chr2,start,end)
        vcfmerge.append(''.join(glob.glob(vcfDir + "*%s*gz"%ref)))
    print("count : %s"%str(count))

    



