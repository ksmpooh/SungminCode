'''
python ../new_format.py /BDATA/smkim/KR_allogenomics/00.oridata_plink/pair_merge/JG.KR.KD.merge_updateID_extractFuncVariant_onlyInpair AR_new_list.withmono.txt

cd /BDATA/smkim/KR_allogenomics/05.AR_2025/with_monomorphic/ped_format
mkdir merge
ls *txt |sed 's/.txt//g' | xargs -I{} -P 10 bash -c "Rscript --vanilla /BDATA/smkim/KR_allogenomics/04.immune_cell_signal/03.data.Pair.matching.R {}.txt merge/{}.pairmerge.txt"
python /BDATA/smkim/KR_allogenomics/05.AR_2025/AR.mergedf.genotypeAAABBB.py *txt

sed 's/.x/_KR/g' AR.20250117.withMono.genotype.txt_ori| sed 's/KR.KD.//g'|sed 's/.y/_KD/g' |sed 's/\_X/\_chr/g' |less -NS
'''

import sys
import os
from collections import defaultdict

def process_file(file_path):
    """
    파일을 읽고 KBA_ID.y와 KBA_ID.x는 그대로, 나머지 열 이름은 고유하게 변경.
    """
    base_name = os.path.basename(file_path)  # 파일 이름 추출
    key = base_name.split("_")[0]  # 파일 이름에서 키워드 추출 (예: ULBP1)

    with open(file_path, "r") as file:
        lines = file.readlines()

    header = lines[0].strip().split("\t")
    data = [line.strip().split("\t") for line in lines[1:]]

    # 열 이름 변경: X6.150289946.G.A.x -> KEY_X6.150289946.G.A.x
    new_header = [
        f"{key}_{col}" if not col.startswith("KBA_ID") else col
        for col in header
    ]

    return new_header, data

def merge_files(file_list, output_file):
    """
    여러 파일을 병합하고 결과를 저장.
    """
    merged_data = defaultdict(list)
    final_header = ["KBA_ID.y", "KBA_ID.x"]  # 초기 헤더

    for file_path in file_list:
        header, data = process_file(file_path)

        # 파일별 데이터 병합
        for row in data:
            key = (row[0], row[1])  # (KBA_ID.y, KBA_ID.x)를 기준으로 병합
            merged_data[key].extend(row[2:])  # 나머지 열 추가

        # 헤더 병합 (중복 방지)
        final_header.extend([col for col in header[2:] if col not in final_header])

    # 결과 파일 저장
    with open(output_file, "w") as out:
        # 헤더 작성
        out.write("\t".join(final_header) + "\n")

        # 데이터 작성
        for key, values in merged_data.items():
            # 기본 키 값과 데이터를 합쳐 작성
            full_row = list(key) + values
            out.write("\t".join(full_row) + "\n")

    print(f"결과가 {output_file}에 저장되었습니다.")

if __name__ == "__main__":
    input_files = sys.argv[1:]  # 입력 파일 리스트
    output_file = "test.txt"  # 결과 파일 이름

    merge_files(input_files, output_file)