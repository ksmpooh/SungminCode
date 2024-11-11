# 우선순위 목록 (높은 순위에서 낮은 순위로)
consequence_priority = [
    "stop_lost",
    "stop_gained",
    "splicing_acceptor/donor",
    "missense",
    "coding",
    "frameshift",
    "3_prime_UTR",
    "5_prime_UTR",
    "downstream",
    "upstream",
    "gene_fusion",
    "deletion",
    "insertion",
    "intronic",
    "splice_region",
    "start_lost",
    "start_retained",
    "stop_retained",
    "synonymous",
    "other"
]

# VEP 파일을 라인 단위로 읽어들임
file_path = 'catalog.GRCh38.with_adjacent_repeats.TRGT.vep'  # 파일 경로를 지정하세요
with open(file_path, 'r') as file:
    lines = file.readlines()

# 결과를 저장할 리스트
result_lines = []

# 데이터 처리
for line in lines:
    columns = line.strip().split('\t')
    
    # "Uploaded_variation"과 "Consequence" 열의 인덱스를 찾음
    if lines.index(line) == 0:
        uploaded_variation_idx = columns.index("Uploaded_variation")
        consequence_idx = columns.index("Consequence")
        result_lines.append(line.strip() + '\tSelected_Consequence\n')
        continue
    
    uploaded_variation = columns[uploaded_variation_idx]
    consequences = columns[consequence_idx].split(',')
    
    # 우선순위에 따라 Consequence를 선택
    selected_consequence = None
    for priority in consequence_priority:
        if priority in consequences:
            selected_consequence = priority
            break
    
    if selected_consequence is None:
        selected_consequence = "other"
    
    # 결과에 추가
    result_line = line.strip() + '\t' + selected_consequence + '\n'
    result_lines.append(result_line)

# 결과를 새로운 파일에 저장
output_file_path = 'selected_consequences_output.txt'  # 원하는 출력 파일 경로를 지정하세요
with open(output_file_path, 'w') as output_file:
    output_file.writelines(result_lines)

print(f"처리가 완료되었습니다. 결과는 {output_file_path}에 저장되었습니다.")