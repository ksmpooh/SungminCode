######
#tool -g gen.gz -snp-stats -osnp output -threads 4
import glob,os

Dir = "/ADATA/smkim/JG/Imputation/"
inDir = "05.imputation/INPUTs/"
outDir = "05.imputation/OUTPUTs/"
shDir= "05.imputation/SCRIPTs/04.gen2vcf/"

mergeDir = outDir + "02.mergeGen/"
refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

gen2vcf = "gen2vcf"
#JG.phasing.chr9.5001_10000.sampleID.sample
def mk_sh(gen):
	print("gen file = " + gen)
	gen_input = mergeDir + gen
    vcf_output = outDir+"04.gen2vcf/" + gen.replace("gen.gz","vcf.gz")
    sample = mergeDir + gen.replace("gen.gz","sample")

    chr = gen.replace("JG.imputation.mergeGen.","").replace("chr","").replace(".gen.gz","").split(".")[0]

	with open(shDir + gen.replace("gen.gz","") + ".gen2vcf.sh",'w') as shwrite:
        shwrite.write("gen2vcf --gen-file %s --gz --sample-file %s --chr %s --out %s"%(gen_input,sample,chr,vcf_output)) 


def main():
        print("main : ...")

	gens = glob.glob(mergeDir+"*processing*gz")
	gen_index = [s.replace(mergeDir,"") for s in gens]

	for gen in gen_index:
		mk_sh(gen)

main()

