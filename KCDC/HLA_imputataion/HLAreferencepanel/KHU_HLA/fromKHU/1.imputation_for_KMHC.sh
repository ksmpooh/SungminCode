#!/bin/bash

# 코드 설명 : input data phasing하는 툴은 Eagle v2.4.1 (가장 최신), imputation 툴은 impute5 v1.2.0 (가장 최신) 을 설치하고 써주시면 됩니다..!

# input genotyp data phasing
inputpath="~~~"
outpath="~~~"

/kimlab_wd/yuo1996/tools/Eagle_v2.4.1/eagle --geneticMapFile /kimlab_wd/yuo1996/tools/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz  \
                                            --numThreads 200 \
                                            --vcf ${inputpath}/input.vcf.gz \
                                            --outPrefix ${outpath}/chr6.MHC24to36mb.eagle.phased
zcat ${outpath}/chr6.MHC24to36mb.eagle.phased.vcf.gz | bgzip > ${outpath}/chr6.MHC24to36mb.eagle.phased.vcf.bgzip.gz
tabix -p vcf ${outpath}/chr6.MHC24to36mb.eagle.phased.vcf.bgzip.gz


# AC field add for IMPUTE5 running
bcftools +fill-tags ${outpath}/chr6.MHC24to36mb.eagle.phased.vcf.bgzip.gz -Oz -o ${outpath}/chr6.MHC24to36mb.eagle.phased.vcf.AC.bgzip.gz -- -t AN,AC
tabix -p vcf ${outpath}/chr6.MHC24to36mb.eagle.phased.vcf.AC.bgzip.gz


refpath="~~" # reference 있는 path
imputepath="~~" # imputation 결과 path

# imputation
mkdir $imputepath/wholesample_imputation_result


/kimlab_wd/yuo1996/tools/impute5/impute5_v1.2.0/impute5_v1.2.0_static \
--h ${refpath}/whole.eagle.phased.sampleQC.AC.xcf.vcf.gz.bcf \
--g ${outpath}/chr6.MHC24to36mb.eagle.phased.vcf.AC.bgzip.gz \
--r 6:24000000-36000000 \
--buffer-region 6:24000000-36000000 \
--m /kimlab_wd/yuo1996/tools/shapeit4/maps/chr6.b38.gmap.gz \    # 첨부드렸습니다, 주석제거 후 돌려주세요
--threads 200 \
--o ${imputepath}/wholesample_imputation_result/wholesample.imputed.vcf.gz \
--l ${imputepath}/wholesample_imputation_result/imputation.log