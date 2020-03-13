######
#tool -g gen.gz -snp-stats -osnp output -threads 4
import glob,os

#Dir = "/DATA/smkim/JG/05.imputation/"
Dir = "/ADATA/smkim/JG/05.imputation/"
inDir = Dir + "INPUTs/"
outDir = Dir + "OUTPUTs/"
shDir = Dir + "SCRIPTs/03.infoscore/"

phasingDir = "/DATA/smkim/Gastric/Phasing/OUTPUTs/03.5Ksplit/"
mergeDir = outDir + "02.mergeGen/"
refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

tool ="/ADATA/smkim/JG/TOOLs/qctool"

#Gastric.phasing.chr9.5001_10000.sampleID.sample
def mk_sh(gen):
	print("gen file = " + gen)
	gen_input = mergeDir + gen
	info_output = outDir+"03.info/" + gen.replace("gen.gz","") + "info"
	with open(shDir + gen.replace("gen.gz","") + "info.sh",'w') as shwrite:
		shwrite.write("%s -g %s -snp-stats -osnp %s"%(tool,gen_input,info_output))


def main():
        print("main : ...")

	gens = glob.glob(mergeDir+"*pro*gz")
	gen_index = [s.replace(mergeDir,"") for s in gens]

	for gen in gen_index:
		mk_sh(gen)

main()


