#minimac4 --chr chr2 --start 1 --end 20000000 --minRatio 0.000001 
#--window 500000 --refhaps reference_panel.chr2.m3vcf.gz 
#--haps chr2.1-20000000.phased.vcf.gz --noPhoneHome --allTypedSites 
#--format GT, DS, GP --prefix chr2.1-20000000.impute 
#--mapFile chr1_geneticmap.txt --referenceEstimates





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
						shwrite.write("%s --chr %s --start %s --end %s --minRatio 0.000001 --window 1000000 --refhaps %s --haps %s --noPhoneHome --allTypedSites --format GT,DS,GP --prefix %s --mapFile %s --referenceEstimates --cpu 1\n"%(minimac4,chr.replace("chr",""),front,tail,ref,VCFin,VCFout,mapin))
						#--cpu
					break


main()


