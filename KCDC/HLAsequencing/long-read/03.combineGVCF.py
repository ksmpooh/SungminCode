#docker run -v "/DATA/smkim/HLA_seq/long-read/03.variant.calling":"/input" -v "/DATA/smkim/HLA_seq/long-read/04.joint.calling":"/output" -v "/DATA/smkim/HLA_seq/REF":"/ref" \
#broadinstitute/gatk gatk CombineGVCFs -R /ref/HLA.target.fasta --variant /input/gvcf.list.txt -O /output/HLA.Longread.Seq.NIH19KT3813_mapped.Q20.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.gvcf.gz


import os,glob

inDir = "/DATA/smkim/HLA_seq/long-read/03.variant.calling"
outDir = "/DATA/smkim/HLA_seq/long-read/04.joint.calling"
refDir = "/DATA/smkim/HLA_seq/REF"


def CombineGVCFs():
    gvcfs = open(inDir + "/gvcf.list.txt","r")
    #gvcfs = [s.replace(inDir,"/input/") for s in gvcfs]
    gvcfs = [s.replace("\n","") for s in gvcfs]

    cmd = "docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" broadinstitute/gatk gatk CombineGVCFs -R /ref/HLA.target.fasta "%(inDir,outDir,refDir)
    for i in gvcfs:
        cmd = cmd + "--variant %s "%i
    
    cmd = cmd + "-O /output/HLA.Longread.Seq.mapped.Q20.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.gvcf.gz"
    print(cmd)
    os.system(cmd)




CombineGVCFs()
#gatk --java-options "-Xmx4g" GenotypeGVCFs \
   #-R Homo_sapiens_assembly38.fasta \
   #-V input.g.vcf.gz \
   #-O output.vcf.gz


def genotypeGVCFs():
    print("genotypeGVCF")
    cmd = "docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" broadinstitute/gatk gatk GenotypeGVCFs -R /ref/HLA.target.fasta "%(inDir,outDir,refDir)
    cmd = cmd + "-V /output/HLA.Longread.Seq.mapped.Q20.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.gvcf.gz "
    cmd = cmd + "-O /output/HLA.Longread.Seq.mapped.Q20.GATK_haplotypeCaller_VariantCalling.GATK_CombineGVCF.genotypeGVCF.vcf.gz"
    print(cmd)
    os.system(cmd)


genotypeGVCFs()

