#shapeit -convert --thread 24 --input-haps Gastric.phasing.chr20 --include-ind id --output-haps test 

import os,glob

outDir = "/ADATA/smkim/JG/04.phasing/OUTPUTs/"
inDir = "/ADATA/smkim/JG/04.phasing/INPUTs/5KsplitSample/"

#a = [12,13]
def shapeit_split(samples):
#	for chr in a:
	for chr in range(5,5+1):
		chr = str(chr)
		print("chrom : "+chr)
		os.system("cp %s02.chr_phasing/JG.phasing.chr%s.haps.gz %s03.5Ksplit/JG.phasing.chr%s.haps.gz"%(outDir,chr,outDir,chr))
		os.system("cp %s02.chr_phasing/JG.phasing.chr%s.sample %s03.5Ksplit/JG.phasing.chr%s.sample"%(outDir,chr,outDir,chr))
		os.system("gzip -d %s03.5Ksplit/JG.phasing.chr%s.haps.gz"%(outDir,chr))
		#os.system("rm %s03.5Ksplit/Gastric.phasing.chr%s.haps"%(outDir,chr))

def main():
	sampleList = glob.glob(inDir+"*.sampleID")
	shapeit_split(sampleList)
	print("Done")
main()

