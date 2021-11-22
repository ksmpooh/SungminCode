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
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/apt-genotype-axiom \
--analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis \
--arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
--dual-channel-normalization true --cel-files $1  --summaries --write-models --out-dir $2/

##SNPolisher
~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-metrics --posterior-file $2/AxiomGT1.snp-posteriors.txt --call-file $2/AxiomGT1.calls.txt --metrics-file $2/AxiomGT1.out.txt

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file $2/AxiomGT1.out.txt --output-dir $2/

~/Downloads/apt-1.19.0-x86_64-intel-linux/bin/ps-classification-supplemental \
--performance-file $2/Ps.performance.txt --summary-file $2/AxiomGT1.summary.txt \
--call-file $2/AxiomGT1.calls.txt --posterior-file $2/AxiomGT1.snp-posteriors.txt \
--output-dir $2/


