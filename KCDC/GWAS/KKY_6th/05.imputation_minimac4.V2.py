import os,glob


minimac4 = "/BDATA/smkim/TOOLs/minimac4"
wdir = "/BDATA/smkim/KKY_6th/"
phasingDir = wdir + "04.phasing/OUTPUTs/02.chr_phasing/"
#phasingDir = "/ADATA/smkim/JG/04.phasing/OUTPUTs/03.5Ksplit/"

imputationDir = wdir + "05.imputation/"


inDir = wdir + "INPUTs/"
inDir = "/BDATA/smkim/JG/05.imputation/INPUTs/"
outDir = imputationDir + "OUTPUTs/01.imputation/"

os.system("mkdir %s"%outDir)
shDir = imputationDir + "SCRIPTs/01.imputationsh/"
os.system("mkdir %s"%shDir)

tool = "/BDATA/smkim/JG/TOOLs/impute4.1.2_r300.2"


# chr21_39000001_44000000.m3vcf.gz
# genetic_map_chr22_combined_b37_addCHR.txt
#refDir = inDir + "KGP3KRG/"
refDir = "/BDATA/smkim/REF/KRG1KGP/m3vcf/"
mapDir = "/BDATA/smkim/REF/map/m3vcf/"


# KKY.6th.phasing.chr16.vcf.gz
def main():
	print("main : ...")

	refs = glob.glob(refDir + "*.gz")
	inputs = glob.glob(phasingDir + "*.vcf.gz")
	count = 0
	for ref in refs:
		count = count + 1
		print("Chunk count : %s"%(str(count)))
		VCFin = ""
		chr,front,tail = ref.replace(refDir,"").replace(".m3vcf.gz","").split("_")
		for fileIn in inputs:
			if VCFin != "":
				break
			tmps = fileIn.replace(phasingDir,"").split(".")
			for tmp in tmps:
				if tmp == chr:
					VCFin = fileIn
					VCFout = VCFin.replace(phasingDir,outDir).replace(".vcf.gz",".%s_%s"%(front,tail)).replace("phasing","imputation_MINIMAC4")
					mapin = mapDir + "genetic_map_%s_combined_b37_addCHR.m3vcf.txt"%chr
					with open(shDir + "Phased.to.imputation.%s.%s_%s.minimca4.sh"%(chr,front,tail),"w") as shwrite:
						shwrite.write("%s --ignoreDuplicates --chr %s --start %s --end %s --minRatio 0.000001 --window 1000000 --refhaps %s --haps %s --noPhoneHome --allTypedSites --format GT,DS,GP --prefix %s --mapFile %s --referenceEstimates --cpu 1\n"%(minimac4,chr.replace("chr",""),front,tail,ref,VCFin,VCFout,mapin))
					break


chunk_ref = "/BDATA/smkim/JG/05.imputation/INPUTs/imputation.POS.auto_forMINIMAC_20220320.txt"
#main()

def main2():
    print("main : ...")

    chunk_list = open(chunk_ref,"r")
    chunks = [s.replace("\n","") for s in chunk_list]
    chunks.pop(0)
    for chunk in chunks:
        ref,imp = chunk.split("\t")
        ref_panel = refDir + "%s.m3vcf.gz"%(ref)
        ref_chr,ref_start,ref_end = ref.split("_")
        imp_chr,imp_start,imp_end = imp.split("_")
        VCFin = glob.glob(phasingDir + "*%s.*vcf.gz"%ref_chr)
        mapin = mapDir + "genetic_map_%s_combined_b37_addCHR.m3vcf.txt"%imp_chr
        #if len(VCFin)
        VCFin = VCFin.pop()
        VCFout = VCFin.replace(phasingDir,outDir).replace(".vcf.gz",".%s_%s"%(imp_start,imp_end)).replace("phasing","imputation_MINIMAC4")
        with open(shDir + "Phased.to.imputation.%s.%s_%s.minimca4.sh"%(imp_chr,imp_start,imp_end),"w") as shwrite:
            shwrite.write("%s --ignoreDuplicates --chr %s --start %s --end %s --minRatio 0.000001 --window 1000000 --refhaps %s --haps %s --noPhoneHome --allTypedSites --format GT,DS,GP --prefix %s --mapFile %s --referenceEstimates --cpu 1\n"%(minimac4,imp_chr.replace("chr",""),imp_start,imp_end,ref_panel,VCFin,VCFout,mapin))


main2()





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

check()

