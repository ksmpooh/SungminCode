#!/usr/bin/env bash
set -euo pipefail

DB_ROOT="/BDATA/smkim/TOOLs/flagger/misc"   # 수정

DB_UNION="00.db_union.bed"
: > "${DB_UNION}"

add_db () {  # chr start end label
  local f="$1"; local prefix="$2"
  [[ -s "$f" ]] || return 0
  local base; base="$(basename "$f" .bed)"
  awk -v L="${prefix}:${base}" 'BEGIN{OFS="\t"}{print $1,$2,$3,L}' "$f" >> "${DB_UNION}"
}

# --- censat 계열 ---
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_bsat.bed"   "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_gsat.bed"   "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_hsat1A.bed" "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_hsat1B.bed" "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_hsat2.bed"  "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_hsat3.bed"  "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_hor.bed"    "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_mon.bed"    "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_only_ct.bed" "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_no_ct.bed"   "censat"
add_db "${DB_ROOT}/stratifications/censat/chm13v2.0_censat.bed"  "censat"

# --- RepeatMasker 레이어 ---
add_db "${DB_ROOT}/stratifications/repeat_masker/chm13v2.0_RM_4.1.2p1_le6_STR.bed"  "STR"
add_db "${DB_ROOT}/stratifications/repeat_masker/chm13v2.0_RM_4.1.2p1_ge7_VNTR.bed" "VNTR"

# --- SD: all 만 사용 ---
add_db "${DB_ROOT}/stratifications/sd/chm13v2.0_SD.all.bed" "sd"

# 정렬 + 완전 동일행 dedup
LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k4,4 "${DB_UNION}" | awk '!seen[$0]++' > "${DB_UNION}.uniq"
mv "${DB_UNION}.uniq" "${DB_UNION}"
echo "[OK] Built ${DB_UNION}"



#!/usr/bin/env bash
set -euo pipefail

LIST="list.txt"  # "KPPD001.1 KPPD001" 형식
LIFT_DIR="01.T2TCHM13_liftover_VS_flagger_dupRegion"
PAF_ANNOT_DIR="02.paf_on_ref_annot"     # (선택) 있으면 MAPQ 붙임
DB_UNION="00.db_union.bed"

OUT_TMP="04.tmp"
OUT_DIR="04.annot_final"
mkdir -p "${OUT_TMP}" "${OUT_DIR}"

# 겹침을 '레이블별 합계 bp'로 요약하고, 최다(overlap bp) 라벨/비율 선택
# MIN_FRAC를 0.0~1.0로 설정하면 '그 이상 겹친 라벨들'도 list로 남김(예: 0.9 = 90%)
MIN_FRAC=${MIN_FRAC:-0.0}

process_one () {
  local id="$1"; local sample="$2"

  local dup="${LIFT_DIR}/${id}.T2TCHM13_liftover_VS_flagger_dupRegion.bed"
  [[ -s "$dup" ]] || { echo "[WARN] missing dup bed: $dup"; return 0; }

  # 0) 정렬 보장
  LC_ALL=C sort -k1,1 -k2,2n -o "$dup" "$dup"

  # 1) (선택) MAPQ/tp/dv 붙이기: 존재하면 map, 없으면 pass-through
  local withmap="${OUT_TMP}/${id}.dup.withMAPQ.bed"
  if [[ -s "${PAF_ANNOT_DIR}/${id}.paf.on_ref.annot.bed" ]]; then
    bedtools map -a "$dup" -b "${PAF_ANNOT_DIR}/${id}.paf.on_ref.annot.bed" \
      -c 6,7,8 -o max,distinct,distinct > "$withmap"
  else
    # Dup만 그대로 통과, MAPQ/tp/dv 자리에 NA 채움
    awk 'BEGIN{OFS="\t"}{print $0,"NA","NA","NA"}' "$dup" > "$withmap"
  fi

  # 2) DB 라벨과 겹침 bp 붙이기
  local labeled="${OUT_TMP}/${id}.dup.withDB.tsv"
  bedtools intersect -a "$withmap" -b "$DB_UNION" -wo > "$labeled" || true

  # 3) 레코드별 요약: 가장 많이 겹친(top) 라벨, top_frac, 그리고 MIN_FRAC 이상 라벨 리스트
  # a 파일 컬럼: [1..N] + [MAPQ_max,tp_list,dv_list] → 여기서는 MAPQ=col5, tp=6, dv=7 (dup 원래 컬럼 개수에 따라 다를 수 있음)
  # labeled: a_cols ... | b_chr b_start b_end label | ovlp_bp
  # 최종 출력: chr start end name score strand MAPQ_max tp_list dv_list top_label top_bp top_frac labels_pass(MIN_FRAC)
  awk -v OFS="\t" -v MINF="${MIN_FRAC}" '
  function flush(){
    if(n==0){return}
    # 총 길이
    dup_len = a_end - a_start
    # top 라벨 찾기 + 리스트 생성
    top_label="NA"; top_bp=0; labels_pass="";
    for(l in sum){
      if(sum[l] > top_bp){ top_bp=sum[l]; top_label=l }
    }
    top_frac = (dup_len>0 ? top_bp/dup_len : 0)
    for(l in sum){
      frac = (dup_len>0 ? sum[l]/dup_len : 0)
      if(frac >= MINF){
        if(labels_pass!="") labels_pass=labels_pass";";
        labels_pass = labels_pass l":"sum[l]":"sprintf("%.4f",frac)
      }
    }
    print a_chr, a_start, a_end, a_name, a_score, a_strand, a_mapq, a_tp, a_dv, top_label, top_bp, sprintf("%.6f",top_frac), labels_pass
  }
  BEGIN{
    # 입력 a(withMAPQ.bed) 예상: chr start end name score strand MAPQ tp dv
    # intersect 결과는 a_컬럼 뒤로 b_chr b_start b_end label ovlp_bp 추가
  }
  {
    # a key = chr|start|end|name|score|strand|MAPQ|tp|dv (고유화하려면 name만 써도 무방)
    key = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9
    if(key != prev){
      if(NR>1) flush()
      # a 정보 보존
      a_chr=$1; a_start=$2; a_end=$3; a_name=$4; a_score=$5; a_strand=$6; a_mapq=$7; a_tp=$8; a_dv=$9;
      # reset
      delete sum; n=0
      prev=key
    }
    # 라벨/겹침 bp 집계 (겹침 없으면 이 줄이 아예 없음 → 마지막 flush 필요)
    label=$(NF-1); bp=$NF
    sum[label]+=bp; n++
  }
  END{
    if(NR>0) flush()
  }' "$labeled" > "${OUT_DIR}/${id}.dup.annot.final.bed"

  echo "[OK] ${id} → ${OUT_DIR}/${id}.dup.annot.final.bed"
}

export -f process_one
export LIFT_DIR PAF_ANNOT_DIR DB_UNION OUT_TMP OUT_DIR MIN_FRAC

# 병렬 처리
xargs -a "$LIST" -n 2 -P 30 bash -c 'process_one "$@"' _