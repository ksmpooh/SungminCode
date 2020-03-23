######
#tool -g gen.gz -snp-stats -osnp output -threads 4
import glob,os

Dir = "/ADATA/smkim/JG/Imputation/"
inDir = "/ADATA/smkim/JG/Imputation/INPUTs/"
outDir = "/ADATA/smkim/JG/Imputation/OUTPUTs/"
shDir= "/ADATA/smkim/JG/Imputation/SCRIPTs/03.info_gen2vcf/"

phasingDir = "/DATA/smkim/JG/Phasing/OUTPUTs/03.5Ksplit/"
mergeDir = outDir + "02.mergeGen/"
refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

qctool ="/ADATA/smkim/JG/TOOLs/qctool"
gen2vcf = "gen2vcf"
#JG.phasing.chr9.5001_10000.sampleID.sample
def mk_sh(gen):
	print("gen file = " + gen)
	gen_input = mergeDir + gen
	info_output = outDir+"03.info/" + gen.replace("gen.gz","") + ".info"
        vcf_output = outDir+"04.gen2vcf/" + gen.replace("gen.gz","vcf.gz")
        sample = mergeDir + gen.replace("gen.gz","sample")

        chr = gen.replace("JG.imputation.mergeGen.","").replace("chr","").replace(".gen.gz","").split(".")[0]

	with open(shDir + gen.replace("gen.gz","") + "info_gen2vcf.sh",'w') as shwrite:
		shwrite.write("%s -g %s -snp-stats -osnp %s\n"%(qctool,gen_input,info_output)) 
                shwrite.write("gen2vcf --gen-file %s --gz --sample-file %s --chr %s --out %s"%(gen_input,sample,chr,vcf_output)) 


def main():
        print("main : ...")

	gens = glob.glob(mergeDir+"*processing*gz")
	gen_index = [s.replace(mergeDir,"") for s in gens]

	for gen in gen_index:
		mk_sh(gen)

main()

