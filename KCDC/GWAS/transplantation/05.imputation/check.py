#JG.imputation.mergeGen.processing.chr8.141000001_146000000.vcf.gz
#JG.imputation.mergeGen.processing.chr11.119000001_124000000.gen.gz
#JG.imputation.mergeGen.processing.chr11.119000001_124000000.gen.gz
#JG.imputation.chr11.55001_60000.sampleID.129000001_134000000.gen.gz

import glob,os

Dir = "/ADATA/smkim/JG/05.imputation/"
inDir = Dir + "INPUTs/"
outDir = Dir + "OUTPUTs/01.imputation/"


check_list = []

def check_region(chr,front,tail):
	#print("check_region... imputation file..made or not")
	static_sample = "1_5000.sampleID"
#	check = outDir + "Gastric.imputation.%s.%s.%s_%s.gen.gz"%(chr,static_sample,front,tail)
        check = outDir + "JG.imputation.%s.%s.%s_%s.gen.gz"%(chr,static_sample,front,tail)
        #check = outDir + "JG.imputation.mergeGen.processing.%s.%s_%s.gen.gz"%(chr,front,tail)
        #check = outDir + "JG.imputation.mergeGen.processing.%s.%s_%s.vcf.gz"%(chr,front,tail)


	if glob.glob(check):
		return
	else:
		check_list.append(chr+"."+front+"_"+tail)


	#file_list = glob.glob(outDir +"*.gen.gz")


def main():
	print("main : ...")

	ref_file = open(inDir + "imputation.POS.auto.txt","r")
	ref_list = [x.replace('\n','').split('\t') for x in ref_file]


	for ref_pos,imputation_pos in ref_list[1:]:
		chr,front,tail = imputation_pos.split('_')
		check_region(chr,front,tail)

main()


outfile = open(Dir+"check_Imputation_position_list.txt", "w")
outfile.write("chr.front_tail\n")

for i in check_list:
	outfile.write(i+"\n")
outfile.close()

