## HLA short-read alignment 
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

wDir = "/DATA/smkim/HLA_seq/"
shortread_Dir = wDir + "MHC_Illumina_Rawdata/"

shortread_raw = shortread_Dir + "Fastq/"
shortread_mapping = shortread_Dir + "01.mapping/"

#os.system("mkdir %s"%(shortread_mapping))

shDir = wDir + "SCRIPTs/"
refDir = "/DATA/smkim/HLA_seq/MHC_Illumina_Rawdata/Fastq/test/"
refDir = "/DATA/smkim/HLA_seq/REF/"
#~/bwa-mem2/bwa-mem2 mem -t 4 ./HLA.target.fasta ../733_S56_L002_R2_001.fastq.gz ../733_S56_L002_R1_001.fastq.gz | samtools sort -o test999.sorted.bam
#bwa mem -M -R "@RG\tID:HWI\tSM:[샘플이름]\tPL:ILLUMINA\tLB:[기계]" -t 10 [참조유전체 fasta] [FASTQ.gz 1; 압축형태 이용가능] [FASTQ.gz 2]
#~/bwa-mem2/bwa-mem2 mem -t 16 -R "@RG\tID:HWI\tSM:733\tPL:ILLUMINA\tLB:Novaseq7000" ./HLA.target.fasta ../733_S56_L002_R2_001.fastq.gz ../733_S56_L002_R1_001.fastq.gz >test_align.sam
#samtools sort test_align.sam -o test_align_sorted.bam

def shortread_mapping_sh(ID_table,col_index):
    print("shortread mapping!")
    tool = "~/bwa-mem2/bwa-mem2"
#    NIH_table = {}
    # ID matching for shortread
   # for i in ID_table:
#        i = i.split()
#        NIH_table[i[col_index]] = i[0]
    shDir2 = shDir + "01.shortread_mapping/"
    os.system("mkdir %s"%shDir2)
    for i in ID_table:
        i = i.split()
        #print(i)
        R1 = shortread_raw + i[2] + ".fastq.gz"
        R2 = R1.replace("_R1_","_R2_")
        #out = shortread_mapping + "HLA.Shortread.Seq.%s.align.sorted.bam"%i[0]
        out = "%sHLA.Shortread.Seq.%s.align.sorted.bam"%(shortread_mapping,i[0])
        with open(shDir2 + "%s_%s.sh"%(i[0],i[2]),"w") as shout:
            shout.write("%s mem -t 16 -R \"@RG\\tID:HWI\\tSM:%s\\tPL:ILLUMINA\\tLB:Novaseq7000\" %s/HLA.target.fasta %s %s | samtools sort -o %s"%(tool,i[0],refDir,R1,R2,out))
        
    #print("Mapping :... sample ID %s -> %s"%())



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
        out = df.replace(".bam",".Deepvariant_VariantCalling")
        with open(shpath,"w") as shout:
            shout.write("sudo docker run -v \"%s\":\"/input\" -v \"%s\":\"/output\" -v \"%s\":\"/ref\" \
                google/deepvariant /opt/deepvariant/bin/run_deepvariant --model_type WGS \
                    -ref /ref/HLA.target.fasta \
                        --reads /input/%s --output_gvcf /output/%s.gvcf.gz \
                            --output_vcf /output/%s.vcf.gz --num_shards 16"%(shortread_mapping,shortread_VC_Dir,refDir,df.replace(shortread_mapping,""),out,out))
    




def main():
    refs = open(wDir + "HLAseq.ID.table.txt","r")
    refs = [s.replace("\n","") for s in refs]

    header = refs.pop(0)
    #shortread_mapping_sh(refs,2)
    deepvariant_calling()

main()