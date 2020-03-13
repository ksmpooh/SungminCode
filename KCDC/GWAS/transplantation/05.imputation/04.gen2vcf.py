######
#tool -g gen.gz -snp-stats -osnp output -threads 4
import glob,os

Dir = "/ADATA/smkim/Gastric/Imputation/"
inDir = "/ADATA/smkim/Gastric/Imputation/INPUTs/"
outDir = "/ADATA/smkim/Gastric/Imputation/OUTPUTs/"
shDir= "/ADATA/smkim/Gastric/Imputation/SCRIPTs/04.gen2vcf/"

phasingDir = "/DATA/smkim/Gastric/Phasing/OUTPUTs/"
mergeDir = outDir + "02_1.mergeGen/"
refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

tool ="gen2vcf"

#Gastric.phasing.chr9.5001_10000.sampleID.sample
def mk_sh(gen):
	print("gen file = " + gen)
	gen_input = mergeDir + gen
	vcf_output = outDir+"04.gen2vcf/" + gen.replace("gen.gz","vcf.gz")
	sample = mergeDir + gen.replace("gen.gz","sample")
	#Gastric.imputation.mergeGen.chr4.132000001_137000000.gen.gz
	chr = gen.replace("Gastric.imputation.mergeGen.","").replace("chr","").replace(".gen.gz","").split(".")[0]
	with open(shDir + gen.replace("gen.gz","") + "gen2vcf.sh",'w') as shwrite:
		shwrite.write("%s --gen-file %s --gz --sample-file %s --chr %s --out %s"%(tool,gen_input,sample,chr,vcf_output)) 


def main():
        print("main : ...")

	gens = glob.glob(mergeDir+"*gz")
	gen_index = [s.replace(mergeDir,"") for s in gens]

#	samples = glob.glob(mergeDir + "*.sample")
#	samples_index = [s for s in samples]
	for gen in gen_index:
		mk_sh(gen)

main()


