## python 02.extract.variant.inGenelist.py [plink] [ref]
## python 02.extract.variant.inGenelist.py KR.KD input.txt


import os
import sys

# 입력 파라미터
plink = sys.argv[1]
ref = sys.argv[2]

# 출력 디렉토리 설정
OUT_DIR = "gene_extract/"
PEDMAP_DIR = "ped_format/"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PEDMAP_DIR, exist_ok=True)

def run_plink(plink, gene, variants):
    """
    PLINK 명령어 실행 및 임시 파일 생성.
    """
    tmp_file = "tmp.txt"
    with open(tmp_file, "w") as tmp_out:
        tmp_out.writelines(f"{var}\n" for var in variants)

    # PLINK 명령 실행
    os.system(f"plink --bfile {plink} --extract {tmp_file} --recodeA --out {OUT_DIR}/KR.KD.{gene}")
    os.system(f"plink --bfile {plink} --extract {tmp_file} --recode ped --out {PEDMAP_DIR}/KR.KD.{gene}")

def convert_ped_to_genotype(gene):
    """
    PED 및 MAP 파일을 읽어 새 포맷으로 변환.
    """
    ped_file = f"{PEDMAP_DIR}/KR.KD.{gene}.ped"
    map_file = f"{PEDMAP_DIR}/KR.KD.{gene}.map"
    output_file = ped_file.replace(".ped", "_covertFormat.genotype.txt")

    # MAP 파일에서 SNP 이름 추출
    with open(map_file, "r") as mapf:
        snp_names = [line.split()[1] for line in mapf]

    # PED 파일 읽고 변환
    with open(ped_file, "r") as pedf:
        output_data = []
        for line in pedf:
            parts = line.strip().split()
            sample_info = parts[:6]
            genotype_data = parts[6:]
            merged_genotypes = ["".join(genotype_data[i:i+2]) for i in range(0, len(genotype_data), 2)]
            output_data.append(sample_info + merged_genotypes)

    # 변환된 결과 저장
    with open(output_file, "w") as outf:
        header = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] + snp_names
        outf.write("\t".join(header) + "\n")
        for row in output_data:
            outf.write("\t".join(row) + "\n")

    print(f"파일이 {output_file}에 저장되었습니다.")

def process_gene_variants(ref_file, plink):
    """
    ref 파일을 읽고 gene별 변이를 처리.
    """
    gene_variants = {}  # gene별 변이 저장
    with open(ref_file, "r") as ref:
        for line in ref:
            var, _, gene = line.strip().split(maxsplit=2)
            gene_variants.setdefault(gene, []).append(var)

    # gene별로 처리
    for gene, variants in gene_variants.items():
        run_plink(plink, gene, variants)
        convert_ped_to_genotype(gene)

def main():
    process_gene_variants(ref, plink)

if __name__ == "__main__":
    main()