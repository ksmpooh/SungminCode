import os,glob


minimac4 = "/BDATA/smkim/TOOLs/minimac4"
wdir = "/BDATA/smkim/BD/"
phasingDir = wdir + "04.phasing/OUTPUTs/02.chr_phasing/"
#phasingDir = "/ADATA/smkim/JG/04.phasing/OUTPUTs/03.5Ksplit/"

imputationDir = wdir + "05.imputation/"

outDir = imputationDir + "OUTPUTs/01.imputation/"
os.system("mkdir %s"%outDir)
shDir = imputationDir + "SCRIPTs/01.imputationsh/"
os.system("mkdir %s"%shDir)

vcfDir = imputationDir + "OUTPUTs/01.imputation/"
vcfMergeDir = imputationDir + "OUTPUTs/02.merge/"
os.system("mkdir %s"%vcfMergeDir)
#shDir = imputationDir + "SCRIPTs/02.merge/"
#os.system("mkdir %s"%shDir)

refDir = "/BDATA/smkim/GWAS/ref/KRG1KGP/m3vcf/"
mapDir = "/BDATA/smkim/GWAS/ref/map/m3vcf/"


# KKY.6th.phasing.chr16.vcf.gz
# KKY.6th.imputation_MINIMAC4"
out_prefix = "BD.2024.imputation_MINIMAC4."


chunk_ref = "/BDATA/smkim/GWAS/imputation.POS.auto_forMINIMAC_20220320.txt"
#main()

def minimac_mksh():
        print("main : ...")
        chunk_list = open(chunk_ref,"r")
        chunks = [s.replace("\n","") for s in chunk_list]
        chunks.pop(0)
        for chunk in chunks:
                ref,imp = chunk.split("\t")
                ref_panel = refDir + "%s.m3vcf.gz"%(ref)
                ref_chr,ref_start,ref_end = ref.split("_")
                imp_chr,imp_start,imp_end = imp.split("_")
                #if chr not in ["chr1","chr2","chr3"]:
                #        continue
                VCFin = glob.glob(phasingDir + "*%s.*vcf.gz"%ref_chr)
                mapin = mapDir + "genetic_map_%s_combined_b37_addCHR.m3vcf.txt"%imp_chr
                #if len(VCFin)
                VCFin = VCFin.pop()
                VCFout = VCFin.replace(phasingDir,outDir).replace(".vcf.gz",".%s_%s"%(imp_start,imp_end)).replace("phasing","imputation_MINIMAC4")
                with open(shDir + "Phased.to.imputation.%s.%s_%s.minimca4.sh"%(imp_chr,imp_start,imp_end),"w") as shwrite:
                        shwrite.write("%s --ignoreDuplicates --chr %s --start %s --end %s --minRatio 0.000001 --window 1000000 --refhaps %s --haps %s --noPhoneHome --allTypedSites --format GT,DS,GP --prefix %s --mapFile %s --referenceEstimates --cpu 1\n"%(minimac4,imp_chr.replace("chr",""),imp_start,imp_end,ref_panel,VCFin,VCFout,mapin))
                        shwrite.write("tabix -f -p vcf %s.dose.vcf.gz"%(VCFout))








def vcf_merge_command(vcflist,chrom,shDir):
        out = open(shDir+chrom+".merge.sh","w")
        out.write("bcftools concat")
        outname = vcfMergeDir+out_prefix+chrom+".vcf.gz"
        for vcf in vcflist:
                out.write(" %s"%vcf)
        out.write(" -Oz >%s\n"%(outname))
        out.write("tabix -f -p vcf %s"%outname)
        out.close()


#KKY.6th.imputation_MINIMAC4.chr22.16050076_21000000.dose.vcf.gz
def vcfMERGE():
        shDir = imputationDir + "SCRIPTs/02.merge/"
        os.system("mkdir %s"%shDir)
        refs = open(chunk_ref,"r")
        refin = [x.replace("\n","") for x in refs]
        refin.pop(0)
        chr1 = "chr1"
        count = 0
        vcfmerge = []
        for refs in refin:
                ref,imp = refs.split("\t")
                chr2,start,end = imp.split("_")
                ref = "%s.%s_%s"%(chr2,start,end)
                if ref == "chr22.51000001_51244237":
                        vcfmerge.append(''.join(glob.glob(vcfDir + "*%s*gz"%ref)))
                        vcf_merge_command(vcfmerge,chr1,shDir)
                        print("%s size: %s"%(chr1,str(len(vcfmerge))))
                        count = count + len(vcfmerge)
                        break
                elif chr1 != chr2:
                        vcf_merge_command(vcfmerge,chr1,shDir)
                        print("%s size: %s"%(chr1,str(len(vcfmerge))))
                        count = count + len(vcfmerge)
                        chr1 = chr2
                        vcfmerge=[]
                ref = "%s.%s_%s"%(chr2,start,end)
                vcfmerge.append(''.join(glob.glob(vcfDir + "*%s*gz"%ref)))
        print("count : %s"%str(count))


#main()
#chr21_39000001_44000000.m3vcf.gz
#Phased.to.imputation.chr2.125000001_130000000.minimca4.sh
def check():
        print("main : ...")
        vcfs = glob.glob(outDir + "*gz")
        shs = glob.glob(shDir + "*.sh")
	#KKY.6th.imputation_MINIMAC4.chr4.35000001_40000000.dose.vcf.gz
        os.system("mkdir %s"%(shDir+"end/"))
        count = 0
        for sh in shs:
                tmp = sh.replace(shDir,"").replace("Phased.to.imputation.","").replace(".minimca4.sh","")
                for vcf in vcfs:
                        if tmp in vcf:
                                os.system("mv %s ./01.imputationsh/end/%s"%(sh,sh.replace(shDir,"")))
                                count = count + 1
                                break
        print("count : %s"%(str(count)))



def main():
        minimac_mksh()
        #vcfMERGE()

main()