######
#tool -no_amf_align -buffer 1000 -int 1 12345 -h reference.hap.gz -l reference.legend.gz -m genetic_map -g phasing.haps.gz -o_gz -o chr1_1_12345
import os,glob

inDir = "/DATA/smkim/Gastric/Imputation/INPUTs/"
#outDir = "/RDATA9/smkim/Gastric/Imputation/"
outDir = "/DATA/smkim/Gastric/Imputation/OUTPUTs/"

shDir = "/DATA/smkim/Gastric/Imputation/SCRIPTs/01.imputationsh/"
phasingDir = "/DATA/smkim/Gastric/Phasing/OUTPUTs/03.5Ksplit/"
tool = "/DATA/smkim/Gastric/TOOLs/impute4.1.2_r300.2"

refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

def mk_sh(samples,chr,front,tail):
	print("mk_sh : ....")
	legend =  refDir+"%s_%s_%s.legend.gz"%(chr,front,tail)
	hap = refDir + "%s_%s_%s.hap.gz"%(chr,front,tail) 
	map = mapDir +"genetic_map_%s_combined_b37.txt"%(chr)
	for sample in samples:
		input = phasingDir + "Gastric.phasing.%s.%s.haps.gz"%(chr,sample)
		output = outDir + "Gastric.imputation.%s.%s.%s_%s"%(chr,sample,front,tail)
		with open(shDir+"Phased.to.imputation"+chr+"."+sample+"."+front+"_"+tail+".sh",'w') as shwrite:
			shwrite.write("%s -no_maf_align -buffer 1000 -int %s %s -h %s -l %s -m %s -g %s -o_gz -o %s"%(tool,front,tail,hap,legend,map,input,output))

def main():
	print("main : ...")
	samples = glob.glob("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/*.sampleID")
	sample_index = [s.replace("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/","") for s in samples]
#	sample_index = samples.replace("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/","")

	preposition = glob.glob(refDir+"chr*_*legend*")
	position = [p.replace(refDir,"").replace(".legend.gz","").split('_') for p in preposition]
	for chr,front,tail in position:
		if chr == "chr20":
			mk_sh(sample_index,chr,front,tail)
main()
