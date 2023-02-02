## HLA short-read alignment  after trmmined
# 20220427
# change ID (set @RG)

# 2277_S83_L002_R2_001.fastq.gz
# input
# ID NIHID
#
#KBA_ID	shortread_ID	Shortread_filename_R1	Shortread_filename_R2	longread_ID	Longread_filename	Longread_filePath
#NIH19KT0247	247	247_S99_L002_R1_001	247_S99_L002_R2_001	2020HLAseq001	1st_Cell.bc1001--bc1001	./MHC_PacBio_Rawdata/1st_Cell_CCS/1st_Cell.bc1001--bc1001
#NIH19KT0248	248	248_S100_L002_R1_001	248_S100_L002_R2_001	2020HLAseq002	1st_Cell.bc1002--bc1002	./MHC_PacBio_Rawdata/1st_Cell_CCS/1st_Cell.bc1002--bc1002
#
import os,glob

#wDir = "/DATA/smkim/HLA_seq/"
wDir = "/BDATA/smkim/HLA_seq/"
#shortread_Dir = wDir + "MHC_Illumina_Rawdata/"
shortread_Dir = wDir + "shortread/"
shDir = shortread_Dir + "SCRIPTs/"
refDir = "/BDATA/smkim/HLA_seq/REF/"

#theme = "hg19"
theme = "hg38_HLAregion"
#theme = "hg38_HLAregion_withALT"
if theme == "hg38_HLAregion_withALT":
    ref = refDir + "hg38/hg38.HLA.region.target.withALT.fasta"
elif theme == "hg38_HLAregion":
    ref = refDir + "hg38/hg38.HLA.region.target.fasta"
else:
    ref = refDir + "HLA.target.fasta"



#shortread_raw = shortread_Dir + "Fastq/"
shortread_raw = wDir+  "00.rawDATA/shortread/"
shortread_mapping = shortread_Dir + "01.mapping_%s/"%(theme)
os.system("mkdir %s"%(shortread_mapping))

#~/bwa-mem2/bwa-mem2 mem -t 4 ./HLA.target.fasta ../733_S56_L002_R2_001.fastq.gz ../733_S56_L002_R1_001.fastq.gz | samtools sort -o test999.sorted.bam
#bwa mem -M -R "@RG\tID:HWI\tSM:[샘플이름]\tPL:ILLUMINA\tLB:[기계]" -t 10 [참조유전체 fasta] [FASTQ.gz 1; 압축형태 이용가능] [FASTQ.gz 2]
#~/bwa-mem2/bwa-mem2 mem -t 16 -R "@RG\tID:HWI\tSM:733\tPL:ILLUMINA\tLB:Novaseq7000" ./HLA.target.fasta ../733_S56_L002_R2_001.fastq.gz ../733_S56_L002_R1_001.fastq.gz >test_align.sam
#samtools sort test_align.sam -o test_align_sorted.bam


#00.rawDATA/trimmed/1003_S71_L002_R1_001_paired.fastq.gz
#00.rawDATA/trimmed/1003_S71_L002_R1_001_unpaired.fastq.gz

def shortread_mapping_sh(ID_table,col_index):
    print("shortread mapping!")
    tool = "/BDATA/smkim/TOOLs/bwa-mem2/bwa-mem2"
#    NIH_table = {}
    # ID matching for shortread
   # for i in ID_table:
#        i = i.split()
#        NIH_table[i[col_index]] = i[0]
    shDir2 = shDir + "01.shortread_mapping_%s/"%(theme)
    os.system("mkdir %s"%shDir2)
    for i in ID_table:
        i = i.split()
        #print(i)
        R1 = shortread_raw + i[2] + "_paired.fastq.gz"
        R2 = R1.replace("_R1_","_R2_")
        #out = shortread_mapping + "HLA.Shortread.Seq.%s.align.sorted.bam"%i[0]
        out = "%sHLA.Shortread.Seq.%s.trimmed.%s_align.sorted.bam"%(shortread_mapping,i[0],theme)
        with open(shDir2 + "%s_%s.sh"%(i[0],i[2]),"w") as shout:
            shout.write("%s mem -t 16 -R \"@RG\\tID:HWI\\tSM:%s\\tPL:ILLUMINA\\tLB:Novaseq7000\" %s %s %s | samtools sort -o %s\n"%(tool,i[0],ref,R1,R2,out))
            shout.write("samtools index %s\n"%(out))
            #shout.write("samtools stats %s > %s.stats\n"%(out,out))
            #shout.write("samtools coverage %s > %s.stats\n"%(out,out))
        
    #print("Mapping :... sample ID %s -> %s"%())

'''
sudo docker run -v "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/Fastq/test/":"/input" -v "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/Fastq/test":"/output" -v "/DATA/smkim/pacbio/INPUTs":"/ref" broadinstitute/picard \
java -jar /usr/picard/picard.jar MarkDuplicates \
I=/input/test_align_sorted.bam \
O=/input/test_align_sorted_dedup.bam \
METRICS_FILE=duplicates \
REMOVE_DUPLICATES=True \
CREATE_INDEX=True
'''

def mark_duplicated():
    print("Mark_duplicated")
    dfs = glob.glob(shortread_mapping + "*.sorted.bam")
    shDir2 = shDir + "01-1.shortread_markDup_%s/"%(theme)
    os.system("mkdir %s"%shDir2)
    for i in dfs:
        i = i.replace(shortread_mapping,"")
        out = i.replace(".bam",".dedup.bam")
        with open(shDir2 + out.replace(".bam",".sh"),"w") as shout:
            shout.write("sudo docker run -v \"%s\":\"/input\" -v \"%s\":\"/ref\" broadinstitute/picard java -jar /usr/picard/picard.jar MarkDuplicates I=/input/%s O=/input/%s METRICS_FILE=duplicates REMOVE_DUPLICATES=True CREATE_INDEX=True\n"%(shortread_mapping,refDir,i,out))
            out = shortread_mapping + out
            shout.write("samtools index %s\n"%(out))
            #shout.write("samtools stats %s > %s.stats\n"%(out,out))
            

        






#sudo docker run -v "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/Fastq/test/":"/input" -v "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/Fastq/test":"/output" -v "/DATA/smkim/pacbio/INPUTs":"/ref" google/deepvariant \
#/opt/deepvariant/bin/run_deepvariant --model_type WGS --ref /ref/HLA.target.fasta \
#--reads /input/test_align_sorted.bam --output_gvcf /output/test_align_sorted_deepvariant.gvcf.gz \
#--output_vcf /output/test_align_sorted_deepvariant.vcf.gz --num_shards 20

def deepvariant_calling():
    print("Variant Calling using Deepvariant")
    shortread_VC_Dir = shortread_Dir + "02.variant.call/"
    os.system("mkdir "+shortread_VC_Dir)

    dfs= glob.glob(shortread_mapping + "*bam")
    shDir2 = shDir + "02.variant.call/"

    for df in dfs:
        shpath = df.replace(shortread_mapping,shDir2).replace(".bam",".DV.variant.calling.sh")
        out = df.replace(".bam",".Deepvariant_VariantCalling").replace(shortread_mapping,"")
        with open(shpath,"w") as shout:
            shout.write("sudo docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" \
                google/deepvariant /opt/deepvariant/bin/run_deepvariant --model_type WGS \
                    -ref /ref/HLA.target.fasta \
                        --reads /input/%s --output_gvcf /output/%s.gvcf.gz \
                            --output_vcf /output/%s.vcf.gz --num_shards 16"%(shortread_mapping,shortread_VC_Dir,refDir,df.replace(shortread_mapping,""),out,out))
    

def gatk_calling():
    print("Variant Calling using GATK haplotypcall")
    shortread_VC_Dir = shortread_Dir + "02.variant.call/"
    os.system("mkdir "+shortread_VC_Dir)

    dfs= glob.glob(shortread_mapping + "*.dedup.bam")
    shDir2 = shDir + "02.variant.call_gatk/"
    os.system("mkdir %s"%shDir2)

    for df in dfs:
        shpath = df.replace(shortread_mapping,shDir2).replace(".bam",".GATK_haplotypeCaller.variant.calling.sh")
        out = df.replace(".bam",".GATK_haplotypeCaller_VariantCalling").replace(shortread_mapping,"")
        with open(shpath,"w") as shout:
            shout.write("sudo docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" \
                -it broadinstitute/gatk bash -c 'gatk HaplotypeCaller \
                    -R /ref/HLA.target.fasta \
                        -I /input/%s -O /output/%s.gvcf.gz \
                            -ERC GVCF'"%(shortread_mapping,shortread_VC_Dir,refDir,df.replace(shortread_mapping,""),out))
    
'''
sudo docker run -v "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/Fastq/test":"/input" -it broadinstitute/gatk \
bash -c 'gatk HaplotypeCaller \
-I /input/test_align_sorted.bam \
-R /input/HLA.target.fasta \
-O /input/test_align_sorted.gatk.gvcf.gz \
-ERC GVCF' \
'''

def main():
    refs = open(wDir + "HLAseq.ID.table_v2.txt","r")
    refs = [s.replace("\n","") for s in refs]

    header = refs.pop(0)

    shortread_mapping_sh(refs,2)
    #mark_duplicated()
    #deepvariant_calling()
    #gatk_calling()

main()