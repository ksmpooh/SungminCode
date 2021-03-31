import os,glob

wDir = "/BDATA/smkim/imputation.tool.check/"
inDir = wDir + "INPUTs/"
outDir = wDir + "OUTPUTs/"
shDir = wDir + "SCRIPTs/impute5/"
#shDir = "/BDATA/smkim/imputation.tool.check/SCRIPTs/impute5/"
os.system("mkdir %s"%(shDir + "02.imputation"))
impDir = outDir + "IMPUTE5/01.imputation/"

#chr1_100000001_105000000.vcf.gz
refDir = inDir + "ref/KGP3KRG/vcf/"

#tool = "/BDATA/smkim/TOOLs/impute5_v1.1.3_static"
tool = "/BDATA/smkim/TOOLs/impute5_1.1.4_static"
#/BDATA/smkim/imputation.tool.check/impute5/impute5_v1.1.3_static \
#--h /BDATA/smkim/imputation.tool.check/impute5/INPUTs/ref/chr1.merge.vcf.gz \
#--g /BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf/defualt.thread.DS.phasing.test.chr1_10Ksample.vcf.gz \
#--o /BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/test/test.vcf.gz \
#--m /BDATA/smkim/imputation.tool.check/impute5/INPUTs/genetic_map_chr1_forIMPUTE5.txt \
#--r 1 --b 1000 --threads 10 \
#--ne 10000
#input = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf/defualt.thread.DS.phasing.test.chr1_10Ksample.vcf.gz"
input = inDir + "phasing/shapeit4/shapeit4.phased.DS.10K.vcf.gz"
#input = "/BDATA/smkim/imputation.tool.check/INPUTs/KBA/Phasing.using.eagle.vcftovcf.KBA.DS.chr1_10Ksample_vcf.vcf.gz"
#input = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/00.haptovcf.defaultealgeoption/eagle.defualtoption.DS.10Ksample.eagle.phasing.chr1_plinkto.eaglehap.vcf.gz"
#input = "/BDATA/smkim/imputation.tool.check/INPUTs/KBA/shapeit4/test.vcf.gz"
#map = "/BDATA/smkim/imputation.tool.check/impute5/INPUTs/genetic_map_chr1_forIMPUTE5.txt"
map = inDir + "map/genetic_map_chr1_forIMPUTE5.txt"
def main():
	ref_file = open(inDir + "imputation.POS.auto.txt","r")
	ref_list = [x.replace('\n','').split('\t') for x in ref_file]
	for ref_pos,imputation_pos in ref_list[1:]:
		ref_chr, ref_front, ref_tail = ref_pos.split("_")
		imp_chr, imp_front, imp_tail = imputation_pos.split("_")
		if ref_chr != "chr1":
			break
#	chr1_10178_5000000.vcf.gz
		ref = refDir + ref_pos + ".vcf.gz"
#	refs = glob.glob(refDir + "chr1_*.vcf.gz")

#	for ref in refs:
		chr = imp_chr
		front = imp_front
		tail = imp_tail
#		chr,front,tail = ref.replace(refDir,"").replace(".vcf.gz","").split("_")
		chr = chr.replace("chr","")
		region = "%s:%s-%s"%(chr,front,tail)
#		out = ref.replace(refDir,impDir).replace(".vcf.gz","impute5.imputation.vcf.gz")
#		out = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/Ref.1kgp.vcf_eagle.phasing/chr1_%s_%s_using.1KGP.haptovcf.ref_eagle.phased.haptovcf.impute5.imputation.vcf.gz"%(front,tail)
#		out = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/Ref.KRG1KGP.haptovcf_eagle.vcftovcf.phased/chr1_%s_%s_using.KRG1KGP.haptovcf.ref_eagle.phased.vcftovcf.impute5.imputation.vcf.gz"%(front,tail)
#		out = ref.replace(refDir,"/ADATA/smkim/imputation.tool.impute5/impute5/").replace(".vcf.gz","impute5.imputation.vcf.gz")
#		out = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/Ref.KRG1KGP.haptovcf_eagle.phased.Defaultoption/chr1_%s_%s_using.KRG1KGP.haptovcf.ref_eagle.phased.haptovcf.impute5.imputation.vcf.gz"%(front,tail)
#		out = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/Ref.KRG1KGP.haptovcf_shapeit4.phasing/chr1_%s_%s_using.KRG1KGP.haptovcf.ref_shapeit4.phased.impute5.imputation.vcf.gz"%(front,tail)
#		out = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/Ref.1kgp.vcf_shapeit4.phasing/chr1_%s_%s_using.1KGP.haptovcf.ref_shapeit4.phased.impute5.imputation.vcf.gz"%(front,tail)
		#out = "/BDATA/smkim/imputation.tool.check/impute5/OUTPUTs/Ref.KRG1KGP.haptovcf_eagle.phased/chr1_%s_%s_using.KRG1KGP.haptovcf.ref_eagle4.phased.phased.imputation.vcf.gz"%(front,tail)
        out = outDir + "Ref.KRG1KGP.haptovcf_shapeit.phased/chr1_%s_%s_using.KRG1KGP.haptovcf.ref_shapeit4.phased.imputation.vcf.gz"%(front,tail)
		with open(ref.replace(refDir,shDir+"02.imputation/").replace(".vcf.gz",".sh"),'w') as shwrite:
#			ref = "/BDATA/smkim/imputation.tool.check/impute5/INPUTs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
			shwrite.write("%s --h %s --g %s --o %s --m %s --r %s --b 1000 --l"%(tool,ref,input,out,map,region,out.replace("vcf.gz","log")))

main()