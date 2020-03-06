import os

tool =   "~/Downloads/Eagle_v2.4.1/eagle"
for chr in range(1,22+1):
	chr = str(chr)
	i= "/ADATA/smkim/JG/04.phasing/OUTPUTs/01.split/JG.merge.chr"+chr
	m = "/ADATA/smkim/JG/04.phasing/INPUTs/map/genetic_map_chr"+chr+"_combined_b37_addCHR.txt"
	os.system(tool + " --bfile "+i+" --geneticMapFile "+m+" --chrom "+chr+" --outPrefix /ADATA/smkim/JG/04.phasing/OUTPUTs/02.chr_phasing/JG.phasing.chr"+chr+" --maxMissingPerSnp 0.3 --maxMissingPerIndiv 0.5 --numThreads 25")



