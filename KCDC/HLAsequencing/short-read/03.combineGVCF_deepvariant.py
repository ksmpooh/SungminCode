

#docker run -v "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/02.variant.call/GATK":"/input" -v "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/03.joint.calling":"/output" -v "/DATA/smkim/HLA_seq/REF":"/ref" \
#broadinstitute/gatk gatk CombineGVCFs -R /ref/HLA.target.fasta -V /input/gvcf.list.txt -O /output/HLA.Shortread.Seq.NIH19KT3814.align.sorted.dedup.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.gvcf.gz
#HLA.Shortread.Seq.NIH19KT3814.align.sorted.dedup.GATK_haplotypeCaller_VariantCalling.gvcf.gz



import os,glob

inDir = "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/02.variant.call/DV"
outDir = "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/03.joint.calling"
refDir = "/DATA/smkim/HLA_seq/REF"


def CombineGVCFs():
    gvcfs = open(inDir + "/gvcf.list.txt","r")
    #gvcfs = [s.replace(inDir,"/input/") for s in gvcfs]
    gvcfs = [s.replace("\n","") for s in gvcfs]

    cmd = "docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" broadinstitute/gatk gatk CombineGVCFs -R /ref/HLA.target.fasta "%(inDir,outDir,refDir)
    for i in gvcfs:
        cmd = cmd + "--variant %s "%i
    
    cmd = cmd + "-O /output/HLA.Shortread.Seq.align.sorted.dedup.Deepvariant_VariantCalling.GATK_CombineGVCF.gvcf.gz"
    #print(cmd)
    os.system(cmd)

#HLA.Shortread.Seq.align.sorted.dedup.Deepvariant_VariantCalling.gvcf.gz

CombineGVCFs()



def genotypeGVCFs():
    print("genotypeGVCF")
    cmd = "docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" broadinstitute/gatk gatk GenotypeGVCFs -R /ref/HLA.target.fasta "%(inDir,outDir,refDir)
    cmd = cmd + "-V /output/HLA.Shortread.Seq.align.sorted.dedup.Deepvariant_VariantCalling.GATK_CombineGVCF.gvcf.gz "
    cmd = cmd + "-O /output/HLA.Shortread.Seq.align.sorted.dedup.Deepvariant_VariantCalling.GATK_CombineGVCF.genotypeGVCF.vcf.gz"
    print(cmd)
    os.system(cmd)


genotypeGVCFs()



