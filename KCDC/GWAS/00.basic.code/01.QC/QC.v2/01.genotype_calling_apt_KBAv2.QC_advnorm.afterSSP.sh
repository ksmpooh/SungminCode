#### apt 1.9 in ~/Download/.....
#### sh 01.geno.sh [cel_files.txt] [SSP_output_dir]

if [ $# -ne 2 ]; then
 echo "Usage: $0 [SSP_output_dir] [Advanced_out_dir]"
 echo "ex) sh 01.geno.sh /DATA/smkim/cel_files.txt /DATA/genocall"
 exit -1
else
 echo "SSP_output_dir : $1"
 echo "Advanced_out_dir : $2"
 echo "ok" 
fi

#/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom
#apt: apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom
#ssp: PN611391_BatchSSP_1.0_Build8/APT/x64-linux/simple-ssps
#adv: advnorm-1.1/advnorm.sh
#lib: Axiom_KBAbetaA.zip (HG38), Axiom_KORV1_1.r1.zip(HG19) 

##GenotypeCalling


/BDATA/smkim/TOOL/advnorm-1.1/advnorm.sh \
  --summary-file $1/AxiomGT1.summary.txt \
  --calls-file $1/AxiomGT1.calls.txt \
  --report-file $1/AxiomGT1.report.txt \
  --trustcheck-file $1/AxiomGT1.trustcheck.txt \
  --analysis-files-path ~/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis \
  --snp-specific-param-file ~/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.snp_specific_analysis.txt \
  --snp-priors-file /BDATA/smkim/TOOL/SSP_model/KBAv1.1_SSP_ver1.models \
  --special-snps-file ~/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.specialSNPs \
  --ps2snp-file ~/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.r1.ps2snp_map.ps \
  --output-dir $2

 # AdvNorm restults to plink format
/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-format-result \
 --calls-file $2/AxiomGT1.calls.txt \
 --annotation-file /home/genome/Downloads/apt-v1.1/Axiom_KORV1.1_Analysis/Axiom_KORV1_1.na35.annot.db \
 --export-plink-file $2/plink/Axiom_KBAv1.1_adv \
 --plink-sort-by-chrpos True \
 --@export-plink-multiallele-snps True \
 --log-file $2/plink/Axiom_KBAv1.1_convert.log

