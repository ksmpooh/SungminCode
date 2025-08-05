#!/usr/bin/env python3

import sys
from cyvcf2 import VCF

def calculate_avg_maxpost(vcf_path):
    vcf = VCF(vcf_path)
    print("VariantID\tAvgMaxPost")  # 헤더 출력

    for variant in vcf:
        gps = variant.format('GP')  # GP 필드 추출
        if gps is None:
            continue  # GP 필드 없으면 건너뛰기

        maxpost_sum = 0.0
        valid_samples = 0

        for sample_gp in gps:
            if sample_gp is None or any(g is None for g in sample_gp):
                continue
            max_gp = max(sample_gp)
            maxpost_sum += max_gp
            valid_samples += 1

        if valid_samples > 0:
            avg_maxpost = maxpost_sum / valid_samples
            print(f"{variant.ID}\t{avg_maxpost:.4f}")
        else:
            print(f"{variant.ID}\tNA")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("사용법: python calculate_avg_maxpost.py <input.vcf.gz>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    calculate_avg_maxpost(vcf_file)


# 헤더 한 줄 추출
#head -1 certainty.KOTRY.AR_2025.chr1.dose.idchange.R20.8.AR_KRonly.MAF0.01_HWE1E-6.txt > certainty.merge.txt

# chr1~chr22 본문을 헤더 제외하고 이어붙이기
#for i in {1..22}; do
  #tail -n +2 certainty.KOTRY.AR_2025.chr${i}.dose.idchange.R20.8.AR_KRonly.MAF0.01_HWE1E-6.txt >> certainty.merge.txt
#done