#!/bin/bash
# conda activate basicGenomics

mkdir 03_assembly
cd 03_assembly

ln -s ../sampleList.txt

awk '$1=="LR"' sampleList.txt | while read -r TYPE PREFIX HIFI; do
    mkdir ${PREFIX} && cd ${PREFIX}

    hifiasm -o ${PREFIX} -t 8 ${HIFI} --hg-size 250m
    # You can specify the estimated genome size using: --hg-size 3g
    # If your sample is inbred, then add: -l 0

    awk '/^S/{print ">"$2;print $3}' ${PREFIX}*.hap1.*p_ctg.gfa > ${PREFIX}.h1.fa
    awk '/^S/{print ">"$2;print $3}' ${PREFIX}*.hap2.*p_ctg.gfa > ${PREFIX}.h2.fa
    cd ..
done

assembly-stats */*.fa > assembly-stat.txt

mkdir DONE
cd ..
