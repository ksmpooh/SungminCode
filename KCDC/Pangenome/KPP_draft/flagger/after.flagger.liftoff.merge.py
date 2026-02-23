#seq 1 132 | awk '{printf "KPPD%03d\n", $1}' | xargs -I {} -P 10 bash -c "python run.2.py /CDATA/pangenome/flagger_test/03.result/{}.winnowmap.sorted.corrected.downsample_1.0.hmm_flagger_prediction.bed /ADATA/pangenome/liftoff/00.rawDATA/{}.1.gencode_GRCh38.liftoff.gff3_polished /ADATA/pangenome/liftoff/00.rawDATA/{}.2.gencode_GRCh38.liftoff.gff3_polished {}.liffoff_with_flaggerregion.txt"


#python run.py final_flagger_prediction.bed KPPD019.1.gencode_GRCh38.liftoff.gff3_polished KPPD019.2.gencode_GRCh38.liftoff.gff3_polished output.txt
#python run.py final_flagger_prediction.bed KPPD019.1.gencode_GRCh38.liftoff.gff3_polished KPPD019.2.gencode_GRCh38.liftoff.gff3_polished output.txt
import sys

def read_flagger(flagger_file):
    """Flagger 파일 읽기 (BED 형식), seqid별로 저장"""
    flagger_dict = {}
    with open(flagger_file, 'r') as f:
        for line in f:
            if line.startswith('track') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            seqid = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            ftype = parts[3]
            if seqid not in flagger_dict:
                flagger_dict[seqid] = []
            flagger_dict[seqid].append((start, end, ftype))
    return flagger_dict

def process_liftoff_genes(liftoff_file, flagger_dict):
    """Liftoff GFF3 파일에서 gene feature만 추출하고 겹치는 flagger feature와 비율 계산"""
    output_lines = []
    with open(liftoff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            if parts[2] != 'gene':
                continue  # gene만 처리

            seqid = parts[0]
            gene_start = int(parts[3]) - 1  # GFF3는 1-based, BED는 0-based
            gene_end = int(parts[4])
            gene_len = gene_end - gene_start

            overlap_dict = {}
            if seqid in flagger_dict:
                for fstart, fend, ftype in flagger_dict[seqid]:
                    # 겹치는 부분 계산
                    overlap_start = max(gene_start, fstart)
                    overlap_end = min(gene_end, fend)
                    if overlap_start < overlap_end:
                        overlap_len = overlap_end - overlap_start
                        if ftype in overlap_dict:
                            overlap_dict[ftype] += overlap_len
                        else:
                            overlap_dict[ftype] = overlap_len

            if overlap_dict:
                # 겹친 길이 기준 내림차순 정렬
                sorted_features = sorted(overlap_dict.items(), key=lambda x: -x[1])
                features_str = ",".join([ftype for ftype, _ in sorted_features])
                ratios_str = ",".join([f"{overlap_len/gene_len:.3f}" for _, overlap_len in sorted_features])
                # 두 컬럼을 맨 앞에 삽입
                parts.insert(0, ratios_str)
                parts.insert(0, features_str)
            else:
                parts.insert(0, "")
                parts.insert(0, "")

            output_lines.append('\t'.join(parts))
    return output_lines

def main():
    if len(sys.argv) != 5:
        print("Usage: python run.py [flagger input] [liftoff.h1.input] [liftoff.h2.input] [output]")
        sys.exit(1)

    flagger_file = sys.argv[1]
    liftoff_h1_file = sys.argv[2]
    liftoff_h2_file = sys.argv[3]
    output_file = sys.argv[4]

    flagger_dict = read_flagger(flagger_file)

    h1_genes = process_liftoff_genes(liftoff_h1_file, flagger_dict)
    h2_genes = process_liftoff_genes(liftoff_h2_file, flagger_dict)

    with open(output_file, 'w') as out:
        for line in h1_genes + h2_genes:
            out.write(line + "\n")

if __name__ == "__main__":
    main()