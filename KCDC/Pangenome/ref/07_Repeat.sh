#!/bin/bash
# conda activate basicGenomics
LINEAGE=angiosperms
THREAD=2
mkdir 06_Repeat
cd 06_Repeat

ln -s ../sampleList.txt

awk '$1=="Assembly" && $4=="haplotype1"' sampleList.txt | while read -r TYPE PREFIX ASM HAP; do

mkdir denovo.${PREFIX} && cd denovo.${PREFIX}

BuildDatabase -name ${PREFIX} ${ASM}
RepeatModeler -database ${PREFIX} -pa ${THREAD} -LTRStruct

cd ..

RepeatMasker -species ${LINEAGE} -s -parallel ${THREAD} -xsmall -alignments ${ASM::-6}.h1.fa
RepeatMasker -species ${LINEAGE} -s -parallel ${THREAD} -xsmall -alignments ${ASM::-6}.h2.fa
RepeatMasker -lib denovo.${PREFIX}/${PREFIX}-families.fa -s -parallel ${THREAD} -xsmall -alignments ${ASM::-6}.h1.fa.masked
RepeatMasker -lib denovo.${PREFIX}/${PREFIX}-families.fa -s -parallel ${THREAD} -xsmall -alignments ${ASM::-6}.h2.fa.masked

mkdir ${PREFIX}.strain-specific && mv *.fa.masked.* ${PREFIX}.strain-specific/
mkdir ${PREFIX}.${LINEAGE} && mv *.fa.* ${PREFIX}.${LINEAGE}

RepeatMasker -species ${LINEAGE} -s -parallel ${THREAD} -alignments ${ASM::-6}.h1.fa
RepeatMasker -species angiosperms -s -parallel ${THREAD} -alignments ${ASM::-6}.h2.fa
RepeatMasker -lib denovo.${PREFIX}/${PREFIX}-families.fa -s -parallel ${THREAD} -alignments ${ASM::-6}.h1.fa.masked
RepeatMasker -lib denovo.${PREFIX}/${PREFIX}-families.fa -s -parallel ${THREAD} -alignments ${ASM::-6}.h2.fa.masked

mkdir ${PREFIX}.strain-specific.hardMasked && mv *.fa.masked.* ${PREFIX}.strain-specific.hardMasked/
mkdir ${PREFIX}.${LINEAGE}.hardMasked && mv *.fa.* ${PREFIX}.${LINEAGE}.hardMasked/

cd ..
done
mkdir DONE
