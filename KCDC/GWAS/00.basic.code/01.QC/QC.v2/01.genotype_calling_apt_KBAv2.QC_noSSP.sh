#### apt 1.9 in ~/Download/.....
#### sh 01.ssp.geno.sh [cel_files.txt] [outDir]

if [ $# -ne 2 ]; then
 echo "Usage: $0 [cel_files.txt] [outDir]"
 echo "ex) sh 01.geno.sh /DATA/smkim/cel_files.txt /DATA/genocall"
 exit -1
else
 echo "cel files : $1"
 echo "outDir : $2"
 echo "ok" 
fi

#/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom
#apt: apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom
#ssp: PN611391_BatchSSP_1.0_Build8/APT/x64-linux/simple-ssps
#adv: advnorm-1.1/advnorm.sh
#lib: Axiom_KBAbetaA.zip (HG38), Axiom_KORV1_1.r1.zip(HG19) 

##GenotypeCalling

/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom \
 --analysis-files-path /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis \
 --arg-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2.xml \
 --do-rare-het-adjustment True \
 --summaries \
 --genotyping-node:snp-priors-input-file /BDATA/smkim/TOOL/SSP_model/KBAv1.1_SSP_ver1.models \
 --genotyping-node:snp-posteriors-output true \
 --artifact-reduction-output-trustcheck true \
 --force TRUE \
 --out-dir $2 \
 --cel-files $1 \
 --log-file $2/log/


/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/ps-metrics \
 --summary-file $2/AxiomGT1.summary.txt \
 --posterior-file $2/AxiomGT1.snp-posteriors.txt \
 --call-file $2/AxiomGT1.calls.txt \
 --metrics-file $2/metrics.txt \
 --log-file $2/ps-metrics.log

/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/ps-classification \
 --species-type human \
 --metrics-file $2/metrics.txt --output-dir $2/classification/ \
 --log-file $2/classification/

/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-format-result \
 --calls-file $2/AxiomGT1.calls.txt \
 --annotation-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.na35.annot.db \
 --export-plink-file $2/plink/Axiom_KBAv1.1_SSP \
 --plink-sort-by-chrpos True \
 --@export-plink-multiallele-snps True \
 --log-file $2/plink/Axiom_KBAv1.1_convert.log



