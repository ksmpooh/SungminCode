##
import os





TOOL = "/BDATA/smkim/TOOLs/minimac4"
SCRIPT = "/BDATA/smkim/imputation.tool.check/SCRIPTs/minimac4/"
OUTPUT = "/BDATA/smkim/imputation.tool.check/OUTPUTs/minimac4/"
MAPS = "/BDATA/smkim/imputation.tool.check/INPUTs/map/genetic_map_chr1_combined_b37_addCHR.txt"
HAPS = "/BDATA/smkim/imputation.tool.check/INPUTs/phasing/shapeit4/shapeit4.phased.DS.10K.vcf.gz"
REFHAPS = "/BDATA/smkim/imputation.tool.check/INPUTs/ref/KGP3KRG/m3vcf/"
#region = open("/BDATA/ghyoon/comparison.imputation.tool/input/imputation.POS.auto_final_20200303.txt",'r')
region = open("/BDATA/smkim/imputation.tool.check/INPUTs/imputation.POS.auto.txt",'r')

region.readline()

while True:
	line = region.readline().replace("\n","")
	if not line: break
	ref_pos,imp_pos = line.split()
	rchr,rfront,rtail = ref_pos.split("_")
	ichr,ifront,itail = imp_pos.split("_")



#	col = line.split()
#	rInfo = col[0].split("_")
	#print rInfo
#	CHROM = rInfo[0]
#	START = rInfo[1]
#	END = rInfo[2]
	if rchr == "chr1":
#		os.system("mkdir /BDATA/ghyoon/comparison.imputation.tool/output/imputation_minimac4/"+col[0])
		os.system("mkdir "+OUTPUT+imp_pos)
		rst = open(SCRIPT+imp_pos+".sh",'w')
#		rst.write(TOOL+" --chr 1 --start "+START+" --end "+END+" --window 1000 --format DS,GT,GP --referenceEstimates --mapFile "+MAPS+" --refHaps "+REFHAPS+str(col[0])+".m3vcf.gz --haps "+HAPS+" --prefix "+OUTPUT+str(col[0])+"/KBA.DS.chr1_10Ksample.minimac4."+str(col[0])+" --log --cpu 1\n")
#		rst.write(TOOL+" --chr 1 --start "+ifront+" --end "+itail+" --window 1000 --format DS,GT,GP --referenceEstimates --mapFile "+MAPS+" --refHaps "+REFHAPS+ref_pos+".m3vcf.gz --haps "+HAPS+" --prefix "+OUTPUT+imp_pos+"_shapit4/KBA.DS.chr1_10Ksample.minimac4."+imp_pos+" --log --cpu 1\n")
#		rst.write(TOOL+" --chr 1 --start "+ifront+" --end "+itail+" --format DS,GT,GP --referenceEstimates --mapFile "+MAPS+" --refHaps "+REFHAPS+ref_pos+".m3vcf.gz --haps "+HAPS+" --prefix "+OUTPUT+imp_pos+"/Default.window.option_KBA.DS.chr1_10Ksample.minimac4."+imp_pos+" --log --cpu 1\n")
#		rst.write(TOOL+" --chr 1 --start "+ifront+" --end "+itail+" --window 1000 --format DS,GT,GP --referenceEstimates --mapFile "+MAPS+" --refHaps "+REFHAPS+ref_pos+".m3vcf.gz --haps "+HAPS+" --prefix "+OUTPUT+imp_pos+"/KBA.DS.chr1_10Ksample.minimac4."+imp_pos+" --log --cpu 1\n")
		rst.write(TOOL+" --chr 1 --start "+ifront+" --end "+itail+" --window 1000 --format DS,GT,GP --referenceEstimates --mapFile "+MAPS+" --refHaps "+REFHAPS+ref_pos+".m3vcf.gz --haps "+HAPS+" --prefix "+OUTPUT+imp_pos+"/shapeit4.phased_KBA.DS.chr1_10Ksample.minimac4."+imp_pos+" --log --cpu 1\n")
		#rst.write(TOOL+" --chr 1 --start "+ifront+" --end "+itail+" --format DS,GT,GP --referenceEstimates --mapFile "+MAPS+" --refHaps "+REFHAPS+ref_pos+".m3vcf.gz --haps "+HAPS+" --prefix "+OUTPUT+imp_pos+"/Default.window.option_shapeit4.phased_KBA.DS.chr1_10Ksample.minimac4."+imp_pos+" --log --cpu 1\n")
		
		rst.close()
region.close()