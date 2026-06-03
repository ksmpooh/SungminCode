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

#/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom
#apt: apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom
#ssp: PN611391_BatchSSP_1.0_Build8/APT/x64-linux/simple-ssps
#adv: advnorm-1.1/advnorm.sh
#lib: Axiom_KBAbetaA.zip (HG38), Axiom_KORV1_1.r1.zip(HG19) 

##GenotypeCalling


/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-genotype-axiom \
 --analysis-files-path /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3 \
 --arg-file /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3/Axiom_KBAbetaB_96orMore_Step2.r2.apt-genotype-axiom.mm.SnpSpecificPriors.AxiomGT1.apt2.xml \
 --out-dir $2 \
 --cel-files $1 \
 --log-file $2/log/ \
 --do-rare-het-adjustment True \
 --allele-summaries True \
 --genotyping-node:snp-posteriors-output True \
 --artifact-reduction-output-trustcheck True



/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/ps-metrics \
 --summary-file $2/AxiomGT1.summary.txt \
 --call-file $2/AxiomGT1.calls.txt \
 --posterior-file $2/AxiomGT1.snp-posteriors.txt \
 --multi-posterior-file $2/AxiomGT1.snp-posteriors.multi.txt \
 --metrics-file $2/metrics.txt \
 --log-file $2/ps-metrics.log

/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/ps-classification \
 --species-type human \
 --metrics-file $2/metrics.txt --output-dir $2/classification/ \
 --multi-metrics-file $2/multi-metrics.txt \
 --log-file $2/classification/

#/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-format-result \
# --calls-file $2/AxiomGT1.calls.txt \
# --annotation-file /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3/Axiom_KBAbetaB.na36.r2.a1.annot.db \
# --export-plink-file $2/plink/Axiom_KBAv2.0B \
# --plink-sort-by-chrpos True \
# --@export-plink-multiallele-snps True \
# --log-file $2/plink/Axiom_KBAv2.0B_convert.log



/BDATA/smkim/TOOL/apt_2.11.8_linux_64_x86_binaries/bin/apt-format-result \
 --calls-file $2/AxiomGT1.calls.txt \
 --annotation-file /BDATA/smkim/TOOL/Axiom_KBAbetaB_r2_3/Axiom_KBAbetaB.na36.r2.a1.annot.db \
 --export-vcf-file $2/vcf/Axiom_KBAv2.0B.vcf \
 --log-file $2/vcf/Axiom_KBAv2.0B_convert.log

bgzip –f $2/vcf/Axiom_KBAv2.0B.vcf

zcat $2/vcf/Axiom_KBAv2.0B.vcf.gz| grep -v "#" | cut -f3,4,5 | grep -P 'N|AX-643521752' | cut -f1 > varN.txt

#zcat $2/vcf/Axiom_KBAv2.0B.vcf.gz | grep -v "#" | cut -f3,4,5 | grep N | cut -f1 > varN.txt
#
bcftools filter -e 'ID=@varN.txt' $2/vcf/Axiom_KBAv2.0B.vcf.gz | bgzip -c > $2/vcf/Axiom_KBAv2.0B.varN.vcf.gz
#
bcftools annotate --header-lines /BDATA/smkim/GWAS/00.genocall_input/header.txt $2/vcf/Axiom_KBAv2.0B.varN.vcf.gz | bcftools annotate --rename-chrs /BDATA/smkim/GWAS/00.genocall_input/1st_rename_chr.txt | bgzip -c > $2/vcf/Axiom_KBAv2.0B.varN.rechr.vcf.gz
#
bcftools norm -m -any $2/vcf/Axiom_KBAv2.0B.varN.rechr.vcf.gz | bgzip -c > $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz
#


#cp $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.annoprb.vcf.gz

bcftools norm -f /BDATA/smkim/HLA_seq/REF/hg38/Homo_sapiens_assembly38.fasta -c x $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz | bcftools annotate --set-id '%ID|%CHROM:%POS:%REF:%ALT' | bgzip -c > $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.lftaln.vcf.gz


plink \
 --vcf $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.lftaln.vcf.gz \
 --double-id \
 --allow-extra-chr \
 --make-bed \
 --out $2/vcf/Axiom_KBAv2.0B.FINAL_genocall

 #plink --bfile $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.lftaln.convert --make-bed --out KOTRY.KR.1stQC

#bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT'  $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz | bgzip -c >  $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.annovar.vcf.gz
#bcftools annotate --set-id '%ID|%CHROM:%POS:%REF:%ALT' $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.vcf.gz | bgzip -c > $2/vcf/Axiom_KBAv2.0B.varN.rechr.norm.annoid.vcf.gz
