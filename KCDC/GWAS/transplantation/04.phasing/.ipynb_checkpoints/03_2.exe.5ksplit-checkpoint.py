#shapeit -convert --thread 24 --input-haps Gastric.phasing.chr20 --include-ind id --output-haps test 

import os,glob


outDir = "/ADATA/smkim/JG/04.phasing/OUTPUTs/"
inDir = "/ADATA/smkim/JG/04.phasing/INPUTs/5KsplitSample/"

shDir = "/ADATA/smkim/JG/04.phasing/SCRIPTs/03.split5k/"
a = [12,13]
def shapeit_split(samples):
#	for chr in a:
	for chr in range(5,5+1):
		chr = str(chr)
		print("chrom : "+chr)
		hap = outDir+"03.5Ksplit/JG.phasing.chr"+chr
		print("hap : " + hap)
		for sample in samples:
			print("sample : ")
			print(sample)
			sample_name = sample.replace(inDir,"")
			print("sample_name : ")
			print(sample_name)
			sample_split = sample_name.split(".")
			print("sample_split : ")
			print(sample_split)
			sampleIn = sample_split[0]+"."+sample_split[1]
			print("sampleIn : ")
			print(sampleIn)
			out = outDir + "03.5Ksplit/JG.phasing.chr"+chr+"."+sampleIn
			print("out :" +out)
			with open(shDir+chr+"_"+sampleIn+".sh",'w') as sh:
				sh.write("shapeit -convert --thread 2 --input-haps %s --include-ind %s --output-haps %s --output-log %s"%(hap,sample,out,out))

def main():
#	hapgzlist = glob.glob(outDir+"2.chr_phasing/*.gz")
	sampleList = glob.glob(inDir+"*.sampleID")
	shapeit_split(sampleList)
	print("Done")
main()

