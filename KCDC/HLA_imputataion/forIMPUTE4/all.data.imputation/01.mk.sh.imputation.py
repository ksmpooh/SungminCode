#tool -no_amf_align -buffer 1000 -int 1 12345 -h reference.hap.gz -l reference.legend.gz -m genetic_map -g phasing.haps.gz -o_gz -o chr1_1_12345
import os,glob


wdir = "/BDATA/smkim/JG/06.HLAimputation/"
inDir = wdir + "INPUTs/"
outDir = wdir + "OUTPUTs/"
shDir = wdir + "SCRIPTs/"
os.system("mkdir "+shDir)

phasingDir = "/BDATA/smkim/JG/04.phasing/OUTPUTs/03.5Ksplit/"

tool = "/BDATA/smkim/JG/TOOLs/impute4.1.2_r300.2"
qctool = "/BDATA/smkim/JG/TOOLs/qctool"
#refDir = inDir + "KGP3KRG/"

#mapDir = inDir + "map/"

def mk_sh(samples,chr,front,tail,ifront,itail):
	shDir = wdir + "SCRIPTs/01.imputationsh/"
	print("mk_sh : ....%schr_%s_%s  to  %s_%s"%(chr,front,tail,ifront,itail))
	#legend =  refDir+"%s_%s_%s.legend.gz"%(chr,front,tail)
    legend = inDir + "Han.legend.gz"
	#hap = refDir + "%s_%s_%s.hap.gz"%(chr,front,tail)
    hap = inDir + "Han.hap.gz"
	#map = mapDir +"genetic_map_%s_combined_b37.txt"%(chr)
    map = inDir + "genetic_map_chr6_combined_b37.txt"
	for sample in samples:
		input = phasingDir + "JG.KR.phasing.%s.%s.haps.gz"%(chr,sample)
		#output = outDir + "JG.KR.imputation.%s.%s.%s_%s"%(chr,sample,ifront,itail)
        output = outDir + "01.HLAimputation/JG.KR.HLAimputation.%s.%s.%s_%s"%(chr,sample,ifront,itail)
		with open(shDir+"Phased.to.imputation"+chr+"."+sample+"."+ifront+"_"+itail+".sh",'w') as shwrite:
			shwrite.write("%s -no_maf_align -int %s %s -h %s -l %s -m %s -g %s -o_gz -o %s"%(tool,ifront,itail,hap,legend,map,input,output))


def merge_sh(samples,chr,front,tail,ifront,itail):
	mergeGenDir = outDir + "01.HLAimputation/"
	with open(shDir+"JG.imputed.%s.%s_%s.sampleMerge.sh"%(chr,ifront,itail),'w') as shwrite:
		shwrite.write("%s "%(qctool))
		for sample in samples:
			sample_input = phasingDir+"JG.KR.phasing.%s.%s.sample"%(chr,sample)
            #gen_input = outDir pre_genDir + "JG.KR.imputation.%s.%s.%s_%s.gen.gz"%(chr,sample,ifront,itail)
			gen_input = outDir + "01.HLAimputation/JG.KR.HLAimputation.%s.%s.%s_%s.gen.gz"%(chr,sample,ifront,itail)
			shwrite.write(" -g %s -s %s"%(gen_input,sample_input))
		output = mergeGenDir + "JG.KR.HLAimputation.mergeGen.%s.%s_%s.gen.gz"%(chr,ifront,itail)
		outsample = mergeGenDir + "JG.KR.imputation.mergeGen.%s.%s_%s.sample"%(chr,ifront,itail)
		shwrite.write(" -og %s -os %s\n"%(output,outsample))
		gen = output
		out = output.replace(".mergeGen.",".mergeGen.processing.")
		shwrite.write("zcat %s |cut -d' ' -f 2- |gzip  -c > %s"%(gen,out))




def main():
	print("main : ...")
	samples = glob.glob("/BDATA/smkim/JG/04.phasing/INPUTs/5KsplitSample/*.sampleID")
	sample_index = [s.replace("/BDATA/smkim/JG/04.phasing/INPUTs/5KsplitSample/","") for s in samples]
#	sample_index = samples.replace("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/","")
	#ref_file = open(inDir + "imputation.POS.auto_final_20200303.txt","r")
	#ref_list = [x.replace('\n','').split('\t') for x in ref_file]
	front = "28477833"
	tail = "33448188"
    #mk_sh(sample_index,"chr6",front,tail,front,tail)
	merge_sh(sample_index,"chr6",front,tail,front,tail)
	

main()