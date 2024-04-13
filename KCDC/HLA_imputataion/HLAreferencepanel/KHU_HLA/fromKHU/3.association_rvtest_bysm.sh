# Rscript ~~~/3_1.association_C4_AA.R --disease ${disease} # 돌리시려는 phenotype 넣어주시면 됩니다. 스크립트 내에서 코드도 수정부탁드립니다.!
#!/bin/bash

disease=$1

#imputed_vcf_path="~~~/wholesample_imputation_result"
outputpath="/BDATA/smkim/JG.HLAimputation/KHU/asso"


# rvtest 설치하시고 돌리시면 됩니다.!
# 아래코드는 제가 예시로 SLE, RA 인 경우 돌아가게 해논 코드 입니다 기본 parameters은 유지해주시고 데이터 phenotype에 맞게 코드 다시 구성하시면 될 것 같아요!
# covariate은 가급적이면 제가 넣어논 것과 동일하게 부탁드립니다! (SEX + genotype PC 1~5)
## QT
mkdir ${outputpath}/$disease
/BDATA/smkim/TOOLs/rvtests/executable/rvtest --inVcf /BDATA/smkim/JG.HLAimputation/KHU/new_20240409/wholesample_imputation_result/wholesample.imputed.vcf.gz \
                                        --out ${outputpath}/$disease/rvtest.$disease \
                                        --covar /BDATA/smkim/JG.HLAimputation/KHU/00.data/pheno/CITY_HLA_QT_pheno_20240408_withPCA_pro.txt \
                                        --covar-name "PC1,PC2,PC3,PC4,PC5" \
                                        --pheno /BDATA/smkim/JG.HLAimputation/KHU/00.data/pheno/CITY_HLA_QT_pheno_20240408_withPCA_pro.txt \
                                        --pheno-name "$disease" \
                                        --dosage "DS" \
                                        --meta "score" \
                                        --freqUpper 0.995000 \
                                        --freqLower 0.005000 \
                                        --numThread 4 \
                                        --noweb



#sed '1d' ${imputed_vcf_path}/DR2.txt | awk '{print $1,$3,$4, $5}'  | tr " " "\t" > ${imputed_vcf_path}/wholesample.imputed.markers
#sed 'ld' §(imputed_vef_path)/DR2.txt l awk '(print $1,$3,$4, $5)' | tr"" "It" > §(imputed_vef_path}/wholesample.imputed.markers
sed '1d' DR2.txt | awk '{print $1,$3,$4, $5}'  | tr " " "\t" > wholesample.imputed.markers
#Rscript ~~~/3_1.association_C4_AA.R --disease ${disease} # 돌리시려는 phenotype 넣어주시면 됩니다. 스크립트 내에서 코드도 수정부탁드립니다.!



#!/bin/bash

disease=$1

#imputed_vcf_path="~~~/wholesample_imputation_result"
outputpath="/SDATA/smkim/KHU/02.asso"


# rvtest 설치하시고 돌리시면 됩니다.!
# 아래코드는 제가 예시로 SLE, RA 인 경우 돌아가게 해논 코드 입니다 기본 parameters은 유지해주시고 데이터 phenotype에 맞게 코드 다시 구성하시면 될 것 같아요!
# covariate은 가급적이면 제가 넣어논 것과 동일하게 부탁드립니다! (SEX + genotype PC 1~5)
## binary
mkdir ${outputpath}/$disease
/BDATA/smkim/TOOLs/rvtests/executable/rvtest --inVcf /SDATA/smkim/KHU/01.imputation/wholesample_imputation_result/wholesample.imputed.vcf.gz \
                                        --out ${outputpath}/$disease/rvtest.$disease \
                                        --covar /SDATA/smkim/KHU/00.data/pheno/CITY_HLA_LOGISTIC_pheno_20240408_withPCA_pro_01to12.txt \
                                        --covar-name "AGE,SEX,PC1,PC2,PC3,PC4,PC5" \
                                        --pheno /SDATA/smkim/KHU/00.data/pheno/CITY_HLA_LOGISTIC_pheno_20240408_withPCA_pro_01to12.txt \
                                        --pheno-name "$disease" \
                                        --dosage "DS" \
                                        --meta "score" \
                                        --freqUpper 0.995000 \
                                        --freqLower 0.005000 \
                                        --numThread 4 \
                                        --noweb


