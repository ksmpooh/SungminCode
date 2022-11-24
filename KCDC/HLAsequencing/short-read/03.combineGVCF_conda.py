import os,glob

inDir = "/BDATA/smkim/HLA_seq/short-read/02.variant.call/gatk"
outDir = "/BDATA/smkim/HLA_seq/short-read/03.joint.calling"
refDir = "/BDATA/smkim/HLA_seq/REF"

#HLA.Shortread.Seq.NIH19KT2304.trimmed.align.sorted.dedup.GATK_haplotypeCaller_VariantCalling.gvcf.gz


def CombineGVCFs():
    gvcfs = open(inDir + "/gvcf.list.txt","r")
    #gvcfs = [s.replace(inDir,"/input/") for s in gvcfs]
    gvcfs = [s.replace("\n","") for s in gvcfs]
    cmd = "gatk CombineGVCFs -R %s "%(refDir + "/HLA.target.fasta")
    #cmd = "docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" broadinstitute/gatk gatk CombineGVCFs -R /ref/HLA.target.fasta "%(inDir,outDir,refDir)
    for i in gvcfs:
        cmd = cmd + "--variant %s/%s "%(inDir,i)
    
    cmd = cmd + "-O %s/HLA.Shortread.Seq.trimmed.align.sorted.dedup.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.gvcf.gz"%(outDir)
    #print(cmd)
    os.system(cmd)


CombineGVCFs()



def genotypeGVCFs():
    print("genotypeGVCF")
    #cmd = "docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" broadinstitute/gatk gatk GenotypeGVCFs -R /ref/HLA.target.fasta "%(inDir,outDir,refDir)
    cmd = "gatk GenotypeGVCFs -R %s "%(refDir + "/HLA.target.fasta")
    cmd = cmd + "-V %s/HLA.Shortread.Seq.trimmed.align.sorted.dedup.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.gvcf.gz "%(outDir)
    cmd = cmd + "-O %s/HLA.Shortread.Seq.trimmed.align.sorted.dedup.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.genotypeGVCF.vcf.gz"%(outDir)
    print(cmd)
    os.system(cmd)


genotypeGVCFs()



