import os

tool =   "~/Downloads/Eagle_v2.4.1/eagle"
for chr in range(1,3+1):
	chr = str(chr)
	i= "/DATA/smkim/Gastric/Phasing/OUTPUTs/splitPlink/gastric.merge.chr"+chr
	m = "/DATA/smkim/Gastric/Phasing/INPUTs/map/genetic_map_chr"+chr+"_combined_b37_addCHR.txt"
	os.system(tool + " --bfile "+i+" --geneticMapFile "+m+" --chrom "+chr+" --outPrefix /DATA/smkim/Gastric/Phasing/OUTPUTs/Gastric.phasing.chr"+chr+" --maxMissingPerSnp 0.3 --maxMissingPerIndiv 0.5 --numThreads 28")



