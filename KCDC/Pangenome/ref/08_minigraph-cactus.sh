#!/bin/bash
# mamba activate basicGenomics
# source /home/junkim/tools/cactus-bin-v2.8.4/venv-cactus-v2.8.4/bin/activate
mkdir 07_Variant && cd 07_Variant

ln -s ../sampleList.txt
awk '$1=="Assembly"{print $2"_"$4"\t"$3}' sampleList.txt >> sequenceFile.tsv

PREFIX=Set1
REF=`head -1 sequenceFile.tsv | awk '{print $1}'`
cactus-pangenome ./jobstorepath ./sequenceFile.tsv --outDir ${PREFIX} --outName ${PREFIX} --reference ${REF} --filter 9 --giraffe clip filter --vcf  --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full --gfa clip full --vcf --giraffe --gfa --gbz --chrom-vg --maxCores 8 --logFile ${PREFIX}.log

mkdir DONE
cd ..


#mkdir -p ${MYBUCKET}/fasta_pp
#cat ${MYBUCKET}/${PREFIX}.seqfile | sed "s\\/data/assembly\\${MYBUCKET}/fasta_pp\\g" > ${MYBUCKET}/${PREFIX}.pp.seqfile
#cactus-preprocess ${MYJOBSTORE} ${MYBUCKET}/${PREFIX}.seqfile ${MYBUCKET}/${PREFIX}.pp.seqfile --maskAlpha --minLength 100000 --brnnCores 16  --realTimeLogging --logFile ${MYBUCKET}/log/${PREFIX}.pp.log
