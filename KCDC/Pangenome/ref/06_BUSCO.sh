#!/bin/bash
#conda activate basicGenomics
mkdir 05_BUSCO
cd 05_BUSCO

for FASTA in `ls -1 ../03_assembly/*/*.fa`; do
    busco -i ${FASTA} -c 8 -o ${FASTA::-3} -m genome --auto-lineage-euk
done
mkdir DONE
