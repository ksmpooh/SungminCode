# one to multi
import os,glob,sys,time

wDir = "/BDATA/smkim/imputation.tool.check/"
inDir = wDir + "INPUTs/"
outDir = wDir + "OUTPUTs/"
shDir = wDir + "SCRIPTs/"
refDir = inDir + "/ref/1KGP/haplegend/"
mapDir = inDir + "map/"


outDir = outDir + "IMPUTE4/V1_5Ksplit_1KGP/"
os.system("mkdir "+outDir)

#python_shDir = outDir + "ex_python/"
#os.system("mkdir "+python_shDir)
#memory_sh = shDir + "check_memory/"
memory_sh = outDir + "check_memory/"
os.system("mkdir "+memory_sh)

logDir = memory_sh + "log/"
os.system("mkdir "+logDir)
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

#	impin1 = inDir + "1_5000.haps.gz"
#	impin2 = inDir + "1_5000.haps.gz"
    impin1 = "/BDATA/smkim/imputation.tool.check/INPUTs/phasing/eagle/split_5K/KBA.DS.chr1.phasing.1_5000.sampleID.haps"
    impin2 = "/BDATA/smkim/imputation.tool.check/INPUTs/phasing/eagle/split_5K/KBA.DS.chr1.phasing.5001_10000.sampleID.haps"
#	sample1 = outDir + "test/test.sample"
#	sample2 = outDir + "test/test.sample"
    sample1 = "/BDATA/smkim/imputation.tool.check/INPUTs/phasing/eagle/split_5K/KBA.DS.chr1.phasing.1_5000.sampleID.sample"
    sample2 = "/BDATA/smkim/imputation.tool.check/INPUTs/phasing/eagle/split_5K/KBA.DS.chr1.phasing.5001_10000.sampleID.sample"

    impute4 = "/BDATA/smkim/JG/TOOLs/impute4.1.2_r300.2"
    qctool ="/BDATA/smkim/JG/TOOLs/qctool"
#	gen2vcf = "/BDATA/smkim/gen2vcf-v2-1.0.4"
    gen2vcf = "gen2vcf-v2-1.0.4"
    for ref_pos,imputation_pos in ref_list[1:]:
        ref_chr, ref_front, ref_tail = ref_pos.split("_")
        imp_chr, imp_front, imp_tail = imputation_pos.split("_")
        if ref_chr != "chr1":
            break
        with open(logDir + "%s.sh"%imputation_pos,"w") as logout:
            logout.write("sh ../memory.check.sh %s_%s"%(imp_front,imp_tail))
        out_chunk = outDir + imputation_pos +"/"
        os.system("mkdir %s"%out_chunk)
	    #legend =  refDir+"%s_%s_%s.legend.gz"%(ref_chr,ref_front,ref_tail)
        legend = refDir + "1000GP_Phase3_chr1.legend.gz"
        hap = refDir + "1000GP_Phase3_chr1.hap.gz"
	    #hap = refDir + "%s_%s_%s.hap.gz"%(ref_chr,ref_front,ref_tail)
        map = mapDir +"genetic_map_%s_combined_b37.txt"%(ref_chr)

#		out = outDir + "test/out/"
        timeout = out_chunk + "time.%s.%s_%s.txt"%(imp_chr,imp_front,imp_tail)
	#	os.system("mkdir %s"%out)
        impout1 = out_chunk + "Tool.imputation.1_5000.%s.%s_%s"%(imp_chr,imp_front,imp_tail)
        impout2 = out_chunk + "Tool.imputation.5001_10000.%s.%s_%s"%(imp_chr,imp_front,imp_tail)
        imp_merge = out_chunk + "Tool.imputation.merge.%s.%s_%s.gen.gz"%(imp_chr,imp_front,imp_tail)
        imp_sample = out_chunk + "Tool.imputation.merge.%s.%s_%s.sample"%(imp_chr,imp_front,imp_tail)
        imp_merge2 = imp_merge.replace(".merge.",".merge.processing.")

        info_output = imp_merge2.replace("gen.gz","info")
        info_output2 = info_output.replace(".info",".forGen2vcf.info")
#		vcfout = imp_merge2.replace("gen.gz","vcf.gz")
        vcfout = imp_merge2.replace(out_chunk,"../%s/"%imputation_pos).replace("gen.gz","vcf.gz")

#		python_shDir = outDir + "ex_python/"
#		os.system("cp %s %s/"%(sample2,out_chunk))
#		with open(shDir + "test.py","w") as shwrite:
        os.system("mkdir "+outDir+"ex.py/")
        with open(outDir + "ex.py/%s.py"%imputation_pos,"w") as shwrite:
            shwrite.write("import os,sys,time\n") # import
            shwrite.write("chunk = \"%s\"\n"%(imputation_pos))
            shwrite.write("a = open(\"%s\",\"w\")\n"%timeout)
            shwrite.write("a.write(\"chunk\ttype\tstart\tend\ttime\\n\")\n")
#imputation sample 1-5000
            shwrite.write("start = time.time()\n") # time start2
            shwrite.write("os.system(\"%s -no_maf_align -buffer 1000 -int %s %s -h %s -l %s -m %s -g %s -o_gz -o %s\")\n"%(impute4,imp_front,imp_tail,hap,legend,map,impin1,impout1)) #
            shwrite.write("end = time.time()\n")
            shwrite.write("shijian = end - start\n")
            shwrite.write("a.write(\"%s\timpute4.1-5000\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
#imputation sample 5001-10000
            shwrite.write("start = time.time()\n") # time start2
            shwrite.write("os.system(\"%s -no_maf_align -buffer 1000 -int %s %s -h %s -l %s -m %s -g %s -o_gz -o %s\")\n"%(impute4,imp_front,imp_tail,hap,legend,map,impin2,impout2)) #
            shwrite.write("end = time.time()\n")
            shwrite.write("shijian = end - start\n")
            shwrite.write("a.write(\"%s\timpute4.5001-10000\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
#merge gen
            shwrite.write("start = time.time()\n") # time start2
            shwrite.write("os.system(\"%s -g %s.gen.gz -s %s -g %s.gen.gz -s %s -og %s -os %s\")\n"%(qctool,impout1,sample1,impout2,sample2,imp_merge,imp_sample))
            shwrite.write("end = time.time()\n")
            shwrite.write("shijian = end - start\n")
            shwrite.write("a.write(\"%s\tmergeGen\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
# gen preprocessing
#zcat merge.gengz | cut -d' ' -f 2- | gzip -c > merge.pro.gen.gz
            shwrite.write("start = time.time()\n") # time start2
            shwrite.write("os.system(\"zcat %s | cut -d ' ' -f 2- | gzip -c > %s\")\n"%(imp_merge,imp_merge2))
            shwrite.write("end = time.time()\n")
            shwrite.write("shijian = end - start\n")
            shwrite.write("a.write(\"%s\tmergeGen.remove.NA\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")


#                        shwrite.write("print(\"time : %s\" + str((time.time() - start)))\n") #
#shwrite.write("%s -g %s -snp-stats -osnp %s\n"%(qctool,gen_input,info_output))
#qctool -g JG.KR.imputation.mergeGen.processing.chr3.128000001_133000000.gen.gz -snp-stats -osnp 000001_133000000.info
# info_score
            shwrite.write("start = time.time()\n")
            shwrite.write("os.system(\"%s -g %s -snp-stats -osnp %s\")\n"%(qctool,imp_merge2,info_output))
            shwrite.write("end = time.time()\n")
            shwrite.write("shijian =  end - start\n")
            shwrite.write("a.write(\"%s\tinfo_score\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
            shwrite.write("os.system(\"grep -v '^#' %s | cut -f1,2,3,4,18 --output-delimiter=' ' > %s\")\n"%(info_output,info_output2))
#gen2vcf
            shwrite.write("start = time.time()\n")
	        #shwrite.write("os.system(\"%s --gen-file %s --info-file %s --sample-file %s --chr %s --out %s\")\n"%(gen2vcf,imp_merge2,info_output2,imp_sample,imp_chr.replace("chr",""),vcfout))
            shwrite.write("os.system(\"./%s --gen-file ../%s --info-file ../%s --sample-file ../%s --chr %s --out %s\")\n"%(gen2vcf,imp_merge2.replace(out_chunk,imputation_pos+"/"),info_output2.replace(out_chunk,imputation_pos+"/"),imp_sample.replace(out_chunk,imputation_pos+"/"),imp_chr.replace("chr",""),vcfout))
                #shwrite.write("%s --g %s --i %s --s %s --c %s --o %s\n"%(tool,gen_input,info_input,sample,chr,vcf_output))
#		shwrite.write("tabix -f -p vcf %s"%vcf_output)
#		shwrite.write("mv %s %s"%(vcf_output,backupDir)
            shwrite.write("end = time.time()\n")
            shwrite.write("shijian =  end - start\n")
            shwrite.write("a.write(\"%s\tgen2vcf\t%s\t%s\t%s\\n\"%(chunk,str(start),str(end),str(shijian)))\n")
            shwrite.write("a.close()")

main()