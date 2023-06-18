import os,glob,sys,time

wDir = "/BDATA/smkim/imputation.tool.check/"
inDir = wDir + "INPUTs/"
outDir = wDir + "OUTPUTs/"
shDir = wDir + "SCRIPTs/"
refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

###ref
#Ref_Pos Imputation_Pos
#chr1_10178_5000000      chr1_10177_5000000
#
###phasing to imputation
#/BDATA/smkim/JG/TOOLs/impute4.1.2_r300.2 -no_maf_align -buffer 1000 -int 99000001 104000000 
#-h /BDATA/smkim/JG/05.imputation/INPUTs/KGP3KRG/chr14_99000001_104000000.hap.gz 
#-l /BDATA/smkim/JG/05.imputation/INPUTs/KGP3KRG/chr14_99000001_104000000.legend.gz 
#-m /BDATA/smkim/JG/05.imputation/INPUTs/map/genetic_map_chr14_combined_b37.txt 
#-g /BDATA/smkim/JG/04.phasing/OUTPUTs/03.5Ksplit/JG.KR.phasing.chr14.75001_77469.sampleID.haps.gz 
#-o_gz -o /BDATA/smkim/JG/05.imputation/OUTPUTs/01.imputation/JG.KR.imputation.chr14.75001_77469.sampleID.99000001_104000000
#shwrite.write("%s -no_maf_align -buffer 1000 -int %s %s -h %s -l %s -m %s -g %s -o_gz -o %s
#shwrite.write("%s -g %s -snp-stats -osnp %s\n"%(qctool,gen_input,info_output)) 

def main():
        ref_file = open(inDir + "imputation.POS.auto.txt","r")
        ref_list = [x.replace('\n','').split('\t') for x in ref_file]

	impin = outDir + "test/test.haps.gz"
	sample = outDir + "test/test.sample"
	impute4 = "/BDATA/smkim/JG/TOOLs/impute4.1.2_r300.2"
	qctool ="/BDATA/smkim/JG/TOOLs/qctool"
	gen2vcf = "/BDATA/smkim/gen2vcf-v2-1.0.4"
	for ref_pos,imputation_pos in ref_list[1:]:
		ref_chr, ref_front, ref_tail = ref_pos.split("_")
		imp_chr, imp_front, imp_tail = imputation_pos.split("_")
		if ref_chr != "chr1":
			break
	        legend =  refDir+"%s_%s_%s.legend.gz"%(ref_chr,ref_front,ref_tail)
	        hap = refDir + "%s_%s_%s.hap.gz"%(ref_chr,ref_front,ref_tail)
	        map = mapDir +"genetic_map_%s_combined_b37.txt"%(ref_chr)

		out = outDir + "test/out/"
		timeout = out + "time.%s.%s_%s.txt"%(imp_chr,imp_front,imp_tail)
	#	os.system("mkdir %s"%out)
		impout = outDir + "test/test.%s.%s_%s"%(imp_chr,imp_front,imp_tail)
		info_output = impout+".info"
		info_output2 = info_output.replace(".info",".forGen2vcf.info")
		vcfout = impout + ".vcf.gz"
		with open(shDir + "test.py","w") as shwrite:
			shwrite.write("import os,sys,time\n") # import
#			shwrite.write("start = time.time()\n") # time start2
			shwrite.write("chunk = \"%s\"\n"%(imputation_pos))
			shwrite.write("a = open(\"%s\",\"w\")\n"%timeout)
                        shwrite.write("a.write(\"chunk\ttype\tstart\tend\ttime\\n\")\n")
#                        shwrite.write("os.system(\"echo hello\")\n") #imputation
                        shwrite.write("start = time.time()\n") # time start2
                        shwrite.write("os.system(\"%s -no_maf_align -buffer 1000 -int %s %s -h %s -l %s -m %s -g %s -o_gz -o %s\")\n"%(impute4,imp_front,imp_tail,hap,legend,map,impin,impout)) #
			shwrite.write("end = time.time()\n")
			shwrite.write("shijian = end - start\n")
			shwrite.write("a.write(\"%s\timpute4\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
#                        shwrite.write("print(\"time : %s\" + str((time.time() - start)))\n") #
#shwrite.write("%s -g %s -snp-stats -osnp %s\n"%(qctool,gen_input,info_output))
#qctool -g JG.KR.imputation.mergeGen.processing.chr3.128000001_133000000.gen.gz -snp-stats -osnp 000001_133000000.info
			shwrite.write("start = time.time()\n")
                        shwrite.write("os.system(\"%s -g %s.gen.gz -snp-stats -osnp %s\")\n"%(qctool,impout,info_output))
			shwrite.write("end = time.time()\n")
                        shwrite.write("shijian =  end - start\n")
                        shwrite.write("a.write(\"%s\tinfo_score\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
			shwrite.write("os.system(\"grep -v '^#' %s | cut -f1,2,3,4,18 --output-delimiter=' ' > %s\")\n"%(info_output,info_output2))
#gen2vcf
                        shwrite.write("start = time.time()\n")
	                shwrite.write("os.system(\"%s --gen-file %s.gen.gz --info-file %s --sample-file %s --chr %s --out %s\")\n"%(gen2vcf,impout,info_output2,sample,imp_chr.replace("chr",""),vcfout))
                #shwrite.write("%s --g %s --i %s --s %s --c %s --o %s\n"%(tool,gen_input,info_input,sample,chr,vcf_output)) 
#		shwrite.write("tabix -f -p vcf %s"%vcf_output)
#		shwrite.write("mv %s %s"%(vcf_output,backupDir)
                        shwrite.write("end = time.time()\n")
                        shwrite.write("shijian =  end - start\n")
                        shwrite.write("a.write(\"%s\tgen2vcf\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
			shwrite.write("a.close()")

main()
