#!/bin/bash
# conda activate basicGenomics

# 1st column = data type
# 2nd column = sample name
# 3rd column = FASTQ files
mkdir 02_LongReadStat
cd 02_LongReadStat

ln -s ../sampleList.txt

echo "sample,length,mean_quality" > ReadLengthQual.csv

awk '$1=="LR"' sampleList.txt | while read -r TYPE PREFIX FASTQ; do
    >${PREFIX}.ReadLengthQual.csv
done
awk '$1=="LR"' sampleList.txt | while read -r TYPE PREFIX FASTQ; do
    pigz -dc -p 100 ${FASTQ} | bioawk -c fastx -v PREFIX="${PREFIX}" 'OFS=","{print PREFIX, length($seq), meanqual($qual)}' >> ${PREFIX}.ReadLengthQual.csv
    # Merge the outputs
    cat ${PREFIX}.ReadLengthQual.csv >> ReadLengthQual.csv
done

cd ..
