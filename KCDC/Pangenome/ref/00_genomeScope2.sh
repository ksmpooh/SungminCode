#!/bin/bash
# conda activate genomescope2
# Create a directory for GenomeScope2 and move into the directory
mkdir 01_genomeSizeEstimation
cd 01_genomeSizeEstimation
ln -s ../sampleList.txt .

# Create a file that contains the path of short-read sequencing data

mkdir tmp

awk '$1=="SR"' sampleList.txt | while read -r TYPE SAMPLE READ1 READ2; do
    echo ${READ1} > pathTo_${SAMPLE}.SRfiles.txt
    echo ${READ2} >> pathTo_${SAMPLE}.SRfiles.txt

    kmc -k21 -t8 -m256 -ci1 -cs10000 @pathTo_${SAMPLE}.SRfiles.txt ${SAMPLE}.reads tmp/
    kmc_tools transform ${SAMPLE}.reads histogram ${SAMPLE}.reads.histo -cx10000
    genomescope2 -i ${SAMPLE}.reads.histo -o ${SAMPLE}.GenomeScope2.output -n ${SAMPLE} -k 21
done

cd ..
