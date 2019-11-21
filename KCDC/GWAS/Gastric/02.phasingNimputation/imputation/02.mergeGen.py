#/DATA/myhwang/TOOLs/qctool_v2.0.1-Ubuntu16.04-x86_64/qctool -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_2_5002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_2_5002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_5002_10002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_5002_10002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_10002_15002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_10002_15002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_15002_20002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_15002_20002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_20002_25002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_20002_25002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_25002_30002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_25002_30002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_30002_35002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_30002_35002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_35002_40002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_35002_40002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_40002_45002.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_40002_45002.sample -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/imputeIMPUTE4/chr1_5000010_6381153_V1_45002_48288.gen.gz -s /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/splitHAPs/chr1_V1_QCed_addINDEL_rmSNP_rmMAF_flipREF_45002_48288.sample -og /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/mergeGEN/chr1_5000010_6381153_V1.gen.gz -os /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/mergeGEN/chr1_5000010_6381153_V1.sample -threads 4
#/DATA/myhwang/TOOLs/qctool_v2.0.1-Ubuntu16.04-x86_64/qctool -g /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/mergeGEN/chr1_5000010_6381153_V1.gen.gz -snp-stats -osnp /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/calINFO/chr1_5000010_6381153_V1_info -threads 4
#-og /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/mergeGEN/chr1_5000010_6381153_V1.gen.gz
#-os /DATA/myhwang/KCHIP_136K/12_checkSNP/RESULTs/mergeGEN/chr1_5000010_6381153_V1.sample -threads 4

import glob,os

Dir = "/DATA/smkim/Gastric/Imputation/"
inDir = "/DATA/smkim/Gastric/Imputation/INPUTs/"
outDir = "/DATA/smkim/Gastric/Imputation/OUTPUTs/"
shDir= "/DATA/smkim/Gastric/Imputation/SCRIPTs/02.mergeGen/"

phasingDir = "/DATA/smkim/Gastric/Phasing/OUTPUTs/03.5Ksplit/"
refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

tool ="/DATA/smkim/Gastric/TOOLs/qctool"

def mk_sh(samples,chr,front,tail):
#	os.system("mkdir "+shDir+chr)
	with open(shDir+chr + "/Gastric.imputed.%s.%s_%s.sampleMerge.sh"%(chr,front,tail),'w') as shwrite:
		shwrite.write("%s "%(tool))
		for sample in samples:
			sample_input = phasingDir+"Gastric.phasing.%s.%s.sample"%(chr,sample)
			gen_input = outDir + "Gastric.imputation.%s.%s.%s_%s.gen.gz"%(chr,sample,front,tail)
			shwrite.write(" -g %s -s %s"%(gen_input,sample_input))
		output = outDir + "mergeGen/Gastric.imputation.mergeGen.%s.%s_%s.hap.gz"%(chr,front,tail)
		outsample = outDir + "mergeGen/Gastric.imputation.mergeGen.%s.%s_%s.sample"%(chr,front,tail)
		shwrite.write(" -og %s -os %s -threads 2"%(output,outsample))
#Gastric.phasing.chr9.5001_10000.sampleID.sample

def main():
        print("main : ...")
        samples = glob.glob("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/*.sampleID")
        sample_index = [s.replace("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/","") for s in samples]
#       sample_index = samples.replace("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/","")

        preposition = glob.glob(refDir+"chr*_*legend*")
        position = [p.replace(refDir,"").replace(".legend.gz","").split('_') for p in preposition]

        for chr,front,tail in position:
                if chr == "chr22":
			print(" chr : %s , front : %s, tail : %s"%(chr,front,tail))
                        mk_sh(sample_index,chr,front,tail)
main()




