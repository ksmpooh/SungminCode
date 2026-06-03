#### apt 1.9 in ~/Download/.....
#### sh 01.geno.sh [SSP_output_dir] [Advanced_out_dir]
adv=$1/../adv
if [ $# -ne 1 ]; then
 echo "Usage: $0 [SSP_output_dir] [Advanced_out_dir]"
 echo "ex) sh 01.geno.sh [SSP_output_dir] [Advanced_out_dir]"
 exit -1
else
 echo "SSP_output_dir : $1"
 echo "Advanced_out_dir : $adv"
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
  --analysis-files-path /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3 \
  --snp-specific-param-file /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3/Axiom_KBAbetaB.r2.probeset_genotyping_parameters.txt \
  --snp-priors-file /BDATA/smkim/TOOL/SSP_model/KBAv2.0B_SSP_v1_20250801.models \
  --special-snps-file /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3/Axiom_KBAbetaB.r2.specialSNPs \
  --ps2snp-file /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3/Axiom_KBAbetaB.r2.ps2snp_map.ps \
  --output-dir $adv

 # AdvNorm restults to plink format
/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-format-result \
 --calls-file $adv/AxiomGT1.calls.txt \
 --annotation-file /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3/Axiom_KBAbetaB.na36.r2.a1.annot.db \
 --export-vcf-file $adv/vcf/Axiom_KBAv2.0B.vcf \
 --log-file $adv/vcf/Axiom_KBAv2.0B_convert.log

bgzip $adv/vcf/Axiom_KBAv2.0B.vcf

zcat $adv/vcf/Axiom_KBAv2.0B.vcf.gz | grep -v "#" | cut -f3,4,5 | grep N | cut -f1 > varN.txt
#
bcftools filter -e 'ID=@varN.txt' $adv/vcf/Axiom_KBAv2.0B.vcf.gz | bgzip -c > $adv/vcf/Axiom_KBAv2.0B.varN.vcf.gz
#
bcftools annotate --header-lines /BDATA/smkim/GWAS/00.genocall_input/header.txt $adv/vcf/Axiom_KBAv2.0B.varN.vcf.gz | bcftools annotate --rename-chrs /BDATA/smkim/GWAS/00.genocall_input/1st_rename_chr.txt | bgzip -c > $adv/vcf/Axiom_KBAv2.0B.varN.rechr.vcf.gz
#
bcftools norm -m -any $adv/vcf/Axiom_KBAv2.0B.varN.rechr.vcf.gz | bgzip -c > $adv/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz
#


#cp $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.annoprb.vcf.gz

bcftools norm -f /BDATA/smkim/HLA_seq/REF/hg38/Homo_sapiens_assembly38.fasta -c x $adv/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz | bcftools annotate --set-id '%ID|%CHROM:%POS:%REF:%ALT' | bgzip -c > $adv/vcf/Axiom_KBAv2.0B.varN.rechr.norm.lftaln.vcf.gz


plink \
 --vcf $adv/vcf/Axiom_KBAv2.0B.varN.rechr.norm.lftaln.vcf.gz \
 --double-id \
 --allow-extra-chr \
 --make-bed \
 --out $adv/vcf/Axiom_KBAv2.0B.FINAL_adv
