import os,sys,time
chunk = "chr1_247000001_249240543"
a = open("/BDATA/smkim/imputation.tool.check/OUTPUTs/test/out/time.chr1.247000001_249240543.txt","w")
a.write("chunk	type	start	end	time\n")
start = time.time()
os.system("/BDATA/smkim/JG/TOOLs/impute4.1.2_r300.2 -no_maf_align -buffer 1000 -int 247000001 249240543 -h /BDATA/smkim/imputation.tool.check/INPUTs/KGP3KRG/chr1_245000001_249240543.hap.gz -l /BDATA/smkim/imputation.tool.check/INPUTs/KGP3KRG/chr1_245000001_249240543.legend.gz -m /BDATA/smkim/imputation.tool.check/INPUTs/map/genetic_map_chr1_combined_b37.txt -g /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.haps.gz -o_gz -o /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543")
end = time.time()
shijian = end - start
a.write("%s	impute4	%s	%s	%s\n"%(chunk,str(start),str(end),str(shijian)))
start = time.time()
os.system("/BDATA/smkim/JG/TOOLs/qctool -g /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543.gen.gz -snp-stats -osnp /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543.info")
end = time.time()
shijian =  end - start
a.write("%s	info_score	%s	%s	%s\n"%(chunk,str(start),str(end),str(shijian)))
os.system("grep -v '^#' /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543.info | cut -f1,2,3,4,18 --output-delimiter=' ' > /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543.forGen2vcf.info")
start = time.time()
os.system("/BDATA/smkim/gen2vcf-v2-1.0.4 --gen-file /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543.gen.gz --info-file /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543.forGen2vcf.info --sample-file /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.sample --chr 1 --out /BDATA/smkim/imputation.tool.check/OUTPUTs/test/test.chr1.247000001_249240543.vcf.gz")
end = time.time()
shijian =  end - start
a.write("%s	gen2vcf	%s	%s	%s\n"%(chunk,str(start),str(end),str(shijian)))
a.close()