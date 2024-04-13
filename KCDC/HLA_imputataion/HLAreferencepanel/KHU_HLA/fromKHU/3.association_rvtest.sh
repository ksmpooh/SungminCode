#!/bin/bash

disease=$1

imputed_vcf_path="~~~/wholesample_imputation_result"
outputpath="~~~"


# rvtest 설치하시고 돌리시면 됩니다.!
# 아래코드는 제가 예시로 SLE, RA 인 경우 돌아가게 해논 코드 입니다 기본 parameters은 유지해주시고 데이터 phenotype에 맞게 코드 다시 구성하시면 될 것 같아요!
# covariate은 가급적이면 제가 넣어논 것과 동일하게 부탁드립니다! (SEX + genotype PC 1~5)

if [[ ${disease} == "SLE" ]]; then

mkdir ${outputpath}/SLE
/home1/rhdfyd/rvtests/executable/rvtest --inVcf ${imputed_vcf_path}/wholesample.imputed.vcf.gz \
                                        --out ${outputpath}/SLE/rvtest.SLE1 \
                                        --covar /kimlab_wd/kah/0.data/1.4th_KCHIP/3.imputed/2.ped/2.4th.kchip.ra.vcfID.ver3.ped \
                                        --covar-name "SEX,PC1,PC2,PC3,PC4,PC5" \
                                        --pheno /kimlab_wd/kah/0.data/1.4th_KCHIP/3.imputed/2.ped/2.4th.kchip.ra.vcfID.ver3.ped \
                                        --pheno-name "SLE_kin" \
                                        --dosage "DS" \
                                        --meta "score" \
                                        --freqUpper 0.995000 \
                                        --freqLower 0.005000 \
                                        --numThread 50 \
                                        --noweb

elif [[ ${disease} == "RA" ]]; then

mkdir ${outputpath}/RA
/home1/rhdfyd/rvtests/executable/rvtest --inVcf ${imputed_vcf_path}/wholesample.imputed.vcf.gz \
                                        --out ${outputpath}/RA/rvtest.RA1 \
                                        --covar /kimlab_wd/kah/0.data/1.4th_KCHIP/3.imputed/2.ped/2.4th.kchip.ra.vcfID.ver3.ped \
                                        --covar-name "SEX,PC1,PC2,PC3,PC4,PC5" \
                                        --pheno /kimlab_wd/kah/0.data/1.4th_KCHIP/3.imputed/2.ped/2.4th.kchip.ra.vcfID.ver3.ped \
                                        --pheno-name "RA_kin" \
                                        --dosage "DS" \
                                        --meta "score" \
                                        --freqUpper 0.995000 \
                                        --freqLower 0.005000 \
                                        --numThread 50 \
                                        --noweb

else
    echo "Invalid value for 'disease'. Please provide either 'SLE' or 'RA'."
fi



sed '1d' ${imputed_vcf_path}/DR2.txt | awk '{print $1,$3,$4, $5}'  | tr " " "\t" > ${imputed_vcf_path}/wholesample.imputed.markers


Rscript ~~~/3_1.association_C4_AA.R --disease ${disease} # 돌리시려는 phenotype 넣어주시면 됩니다. 스크립트 내에서 코드도 수정부탁드립니다.!