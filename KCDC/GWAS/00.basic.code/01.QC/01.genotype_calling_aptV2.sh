#### apt 1.9 in ~/Download/.....
#### sh 01.geno.sh [cel_files.txt] [outDir]

if [ $# -ne 2 ]; then
 echo "Usage: $0 [cel_files.txt] [outDir]"
 echo "ex) sh 01.geno.sh /DATA/smkim/cel_files.txt /DATA/genocall"
 exit -1
else
 echo "cel files : $1"
 echo "outDir : $2"
 echo "ok" 
fi


##GenotypeCalling
~/Downloads/apt_2.11.4_linux_64_bit_x86_binaries/bin/apt-genotype-axiom \
--analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis \
--arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
--do-rare-het-adjustment true \
--summaries -dual-channel-normalization true \
--genotyping-node:snp-posteriors-output true \
--genotyping-node:snp-priors-input-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.generic_prior.txt \
--out-dir $2 \
--cel-files $1 \
--log-file $2/log/

#--summaries -dual-channel-normalization true \
#--genotyping-node:snp-posteriors-output true \



##SNPolisher

~/Downloads/apt_2.11.4_linux_64_bit_x86_binaries/bin/ps-metrics \
--summary-file $2/AxiomGT1.summary.txt \
--posterior-file $2/AxiomGT1.snp-posteriors.txt \
--call-file $2/AxiomGT1.calls.txt \
--metrics-file $2/metrics.txt

~/Downloads/apt_2.11.4_linux_64_bit_x86_binaries/bin/ps-classification --species-type human \
--metrics-file $2/metrics.txt --output-dir $2/classification/ \
--log-file $2/classification/


# conert plink
~/Downloads/apt_2.11.4_linux_64_bit_x86_binaries/bin/apt-format-result \
--calls-file $2/AxiomGT1.calls.txt \
--annotation-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.na35.annot.db \
--export-plink-file $2/plink/Axiom_KBAv1.1 \
--plink-sort-by-chrpos True \
--log-file $2/plink/Axiom_KBAv1.1_convert.log
