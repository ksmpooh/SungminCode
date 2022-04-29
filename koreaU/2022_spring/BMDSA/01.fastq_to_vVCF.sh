#! /bin/bash
#### ./01.fastq_to_VCF.sh [R1.fastq] [R2.fastq]

if [ $# -ne 2 ]; then
 echo "Usage: $0 [R1.fastq or fastq.gz] [R2.fastq or fastq.gz]"
 echo "ex1) ./01.fastq_to_VCF.sh sampleID_xxxx_R1_xxxx_fastq sampleID_xxxx_R2_xxxx_fastq"
 echo "ex2) ./01.fastq_to_VCF.sh fastq/sampleID_xxxx_R1_xxxx_fastq.gz fastq/sampleID_xxxx_R2_xxxx_fastq.gz"
 exit -1
else
 echo "R1 : $1"
 echo "R2 : $2"
 echo "ok" 
fi
#./01.fastq_to_VCF.sh fastq/SNUM0377-01_L001_R1_001.fastq.gz fastq/SNUM0377-01_L001_R2_001.fastq.gz

sampleID=$(echo $1 | tr '/' '\t' | awk '{print $NF}'| cut -d"_" -f1)
echo "Sample ID : $sampleID"


echo "Alignment : Fastq to mapped Map"
bwa.kit/bwa mem -t 8 -M -R "@RG\tID:${sampleID}\tPL:ILLUMINA\tPU:${sampleID}\tSM:$sampleID\tLB:${sampleID}-1" ref/human_g1k_v37_decoy.fasta <(zcat $1) <(zcat $2) \
| /home/ec2-user/programs/samtools-1.10/bin/samtools view -huS \
| /home/ec2-user/programs/samtools-1.10/bin/samtools sort -@ 8 -m 2G -o $sampleID.bam -O bam -T $sampleID.tmp

java -jar /home/ec2-user/bms553_lab08/picard/picard.jar MarkDuplicates I=$sampleID.bam O=dedup.$sampleID.bam M=markdups_$sampleID.txt \
ASSUME_SORT_ORDER=queryname MAX_RECORDS_IN_RAM=2000000 COMPRESSION_LEVEL=1 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
PROGRAM_RECORD_ID=null ADD_PG_TAG_TO_READS=false READ_NAME_REGEX=null; /home/ec2-user/programs/samtools-1.10/samtools index /home/ec2-user/bms553_lab08/dedup.$sampleID.bam


echo "Recalibration"
gatk-4.1.6.0/gatk BaseRecalibrator -R ref/human_g1k_v37_decoy.fasta -I dedup.$sampleID.bam --known-sites ref/dbsnp_138.b37.vcf \
--known-sites ref/1000G_phase1.snps.high_confidence.b37.vcf --known-sites ref/Mills_and_1000G_gold_standard.indels.b37.vcf \
-L ref/hs37d5_MedExome_and_coding_padded_Nov2015.bed -O $sampleID.recal.table

gatk-4.1.6.0/gatk ApplyBQSR -R ref/human_g1k_v37_decoy.fasta -I dedup.$sampleID.bam -O $sampleID.cram -bqsr $sampleID.recal.table \
--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30

/home/ec2-user/programs/samtools-1.10/samtools index /home/ec2-user/bms553_lab08/$sampleID.cram


echo "Variant calling"
gatk-4.1.6.0/gatk HaplotypeCaller -R ref/human_g1k_v37_decoy.fasta -I $sampleID.cram -O $sampleID.g.vcf.gz -ERC GVCF -L ref/hs37d5_MedExome_and_coding_padded_Nov2015.bed
gatk-4.1.6.0/gatk GenotypeGVCFs -R ref/human_g1k_v37_decoy.fasta -V $sampleID.g.vcf.gz -O rawcalls.$sampleID.vcf.gz -L ref/hs37d5_MedExome_and_coding_padded_Nov2015.bed
