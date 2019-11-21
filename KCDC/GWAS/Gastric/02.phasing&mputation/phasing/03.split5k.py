#shapeit -convert --thread 24 --input-haps Gastric.phasing.chr20 --include-ind id --output-haps test 

import os,glob

outDir = "/DATA/smkim/Gastric/Phasing/OUTPUTs/"
inDir = "/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/"

def shapeit_split(samples):
	for chr in range(1,22+1):
		chr = str(chr)
		print("chrom : "+chr)
		os.system("cp %s/02.chr_phasing/Gastric.phasing.chr%s.haps.gz %s/03.5Ksplit/Gastric.phasing.chr%s.haps.gz"%(outDir,chr,outDir,chr))
		os.system("cp %s/02.chr_phasing/Gastric.phasing.chr%s.sample %s/03.5Ksplit/Gastric.phasing.chr%s.sample"%(outDir,chr,outDir,chr))
		os.system("gzip -d %s/03.5Ksplit/Gastric.phasing.chr%s.haps.gz"%(outDir,chr))
#		hap = outDir+"03.5Ksplit/Gastric.phasing.chr"+chr
#		print("hap : " + hap)
#		for sample in samples:
#			print("sample : ")
#			print(sample)
#			sample_name = sample.replace(inDir,"")
#			print("sample_name : ")
#			print(sample_name)
#			sample_split = sample_name.split(".")
#			print("sample_split : ")
#			print(sample_split)
#			sampleIn = sample_split[0]+"."+sample_split[1]
#			print("sampleIn : ")
##			print(sampleIn)
#			out = outDir + "03.5Ksplit/Gastric.phasing.chr"+chr+"."+sampleIn
#			print("out :" +out)
#			os.system("shapeit -convert --thread 28 --input-haps %s --include-ind %s --output-haps %s --output-log %s"%(hap,sample,out,out))
#			os.system("gzip -c %s.haps > %s.haps.gz"%(out,out))
#			os.system("rm %s.haps"%(out))
#		os.system("rm %s/03.5Ksplit/Gastric.phasing.chr%s.haps"%(outDir,chr))

def rm_unziped_hap():
	for chr in range(4,22+1):
		os.system("rm %s/03.5Ksplit/Gastric.phasing.chr%s.haps"%(outDir,chr))

def main():
#	hapgzlist = glob.glob(outDir+"2.chr_phasing/*.gz")
	sampleList = glob.glob(inDir+"*.sampleID")
	shapeit_split(sampleList)
	print("Done")
main()
#rm_unziped_hap()
