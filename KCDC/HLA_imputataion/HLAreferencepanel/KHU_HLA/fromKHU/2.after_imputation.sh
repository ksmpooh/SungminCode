#!/bin/bash

# 코드 설명 : imputation 결과파일 (wholesample.imputed.vcf.gz) 에서 dosage를 추출해낸 후 Diploid C4 copy number dosage를 계산합니다.

# input genotyp data phasing
inputpath="~~~"
outpath="~~~"

cd ~~~/wholesample_imputation_result # imputation 결과 파일 들어있는 폴더.

zcat wholesample.imputed.vcf.gz | grep -v "##" | grep IMP | cut -f 1-2 > imputed_snps.txt
zcat wholesample.imputed.vcf.gz | grep -v "#" | cut -f 1-5 > ref_snps.txt

cd ${imputepath}/wholesample_imputation_result
awk '$3 ~ /^AA|HLA|copy/' ref_snps.txt | cut -f 1-2 > C4HLA_snps.txt

# extract HLA/C4; 좀 걸림.
vcftools --gzvcf wholesample.imputed.vcf.gz \
         --positions C4HLA_snps.txt \
         --recode \
         --recode-INFO-all \
         --stdout | gzip > wholesample.imputed.for_C4HLA_check.vcf.gz

zcat wholesample.imputed.for_C4HLA_check.vcf.gz | java -jar /kimlab_wd/yuo1996/tools/HLA-TAPAS-master/dependency/vcf2beagle.jar 0 wholesample.imputed.for_C4HLA_check
gzip -d wholesample.imputed.for_C4HLA_check.bgl.gz
#gzip -d wholesample.imputed.for_C4HLA_check.vcf.gz
#sed -i '/^##/d' wholesample.imputed.for_C4HLA_check.vcf
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' wholesample.imputed.for_C4HLA_check.vcf.gz -H --output MHC_GT.txt


# diploid C4 copy calculation

cd ~~~/wholesample_imputation_result

bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' wholesample.imputed.vcf.gz -H --output dosage.txt
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%AF\t%INFO\n' wholesample.imputed.vcf.gz -H --output DR2.txt

cat <(head -n 1 dosage.txt) <(cut -f 2 ./imputed_snps.txt | grep -F -f - dosage.txt) > dosage_imputed_snps_only.txt
cat <(head -n 1 DR2.txt) <(cut -f 2 ./imputed_snps.txt | grep -F -f - DR2.txt) > DR2_imputed_snps_only.txt


Rscript ~~~/2_1.mk_diploid_table.R --snpver ${ver} --haplonetver ${haplonetver}

# imputation QC  : DR2 > 0.5 , MAF >  0.005 (0.5%)%