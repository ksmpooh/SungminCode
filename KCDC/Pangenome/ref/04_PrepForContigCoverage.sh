#!/bin/bash
# conda activate basicGenomics
LENGTH=250000000

mkdir 04_NGx
cd 04_NGx

ln -s ../sampleList.txt .

awk '$1=="LR"' sampleList.txt | while read -r TYPE SAMPLE FASTQ; do
    ls -1 ../03_assembly/${SAMPLE}/${SAMPLE}.*.fa | awk -v sample="${SAMPLE}" 'OFS="\t"{print "Assembly",sample,$1,"haplotype"NR}' >> sampleList.txt
done

echo "Sample,Length,Type,Coverage,LengthSum" > length.csv

awk '$1=="Assembly"' sampleList.txt | while read -r TYPE SAMPLE FASTA HAP; do
    cat ${FASTA} | bioawk -c fastx -v sample="${SAMPLE}.${HAP}" 'OFS=","{print sample, length($seq)}' | sort -k2rV -t "," | \
        awk -F "," -v len="${LENGTH}" -v type="${TYPE}" 'OFS=","{ print $1,$2,type,(sum+0)/len,(sum+0); sum+=$2 } END { print $1,$2,type,sum/len,sum }' \
        >> length.csv
done
