#!/usr/bin/env bash
set -euo pipefail

#### usage:
#### bash 02.snpolisher.marker_select.sh [ori_geno_dir] [ssp_geno_dir] [adv_dir] [merge_output_prefix] [KBA_version]

if [ $# -ne 5 ]; then
    echo "Usage: $0 [ori_geno_dir] [ssp_geno_dir] [adv_dir] [merge_output_prefix] [KBA_version]"
    echo "Example: $0 /DATA/genocall /DATA/genocall/SSP /DATA/genocall/adv KOTRY.AR_2025_KRKD.1stQC KBAv1.1"
    echo "Available KBA_version: KBAv1.1 | KBAv2.0A | KBAv2.0B"
    exit 1
fi

ori_geno_dir=$1
ssp_geno_dir=$2
adv_dir=$3
merge_prefix=$4
kba_version=$5

echo "ori_geno_dir        : ${ori_geno_dir}"
echo "ssp_geno_dir        : ${ssp_geno_dir}"
echo "adv_dir             : ${adv_dir}"
echo "merge_output prefix : ${merge_prefix}"
echo "KBA version         : ${kba_version}"
echo

########################################
# functions
########################################

check_file() {
    local f="$1"
    if [ ! -f "$f" ]; then
        echo "[ERROR] File not found: $f" >&2
        exit 1
    fi
}

# probeset_id list -> actual BIM marker ID list
# input:
#   $1 = probeset_id list (AX-xxxx 형식)
#   $2 = bim file
#   $3 = output extract list
#
# bim 2nd col 예:
#   AX-100005102|chr1:69610:C:T
# 여기서 | 앞쪽 AX-100005102 를 probeset_id 와 매칭
make_extract_from_bim() {
    local probeset_list="$1"
    local bim_file="$2"
    local out_list="$3"

    check_file "$probeset_list"
    check_file "$bim_file"

    awk '
        NR==FNR {
            ids[$1]=1
            next
        }
        {
            split($2, a, "|")
            if (a[1] in ids) print $2
        }
    ' "$probeset_list" "$bim_file" > "$out_list"
}

########################################
# 1. SNP selection list 만들기 (classification 기준)
########################################

cat \
  <(cat "${ori_geno_dir}/classification/CallRateBelowThreshold.ps" \
         "${ori_geno_dir}/classification/Other.ps" \
         "${ori_geno_dir}/classification/OffTargetVariant.ps" | grep -v '^probeset_id$') \
  <(cat "${ssp_geno_dir}/classification/Recommended.ps" | grep -v '^probeset_id$') \
  | sort | uniq -c | awk '$1==2{print $2}' \
  > ssp.snp.select.txt

grep -v -F -f \
  <(cat "${ori_geno_dir}/classification/Recommended.ps" \
         "${ssp_geno_dir}/classification/Recommended.ps" \
      | grep -v '^probeset_id$' \
      | sort | uniq) \
  <(cat "${adv_dir}/classification/Recommended.ps" | grep -v '^probeset_id$') \
  > adv.snp.select.txt

grep -v '^probeset_id$' "${ori_geno_dir}/classification/Recommended.ps" > ori.recommended.txt

echo "[INFO] ori.recommended.txt : $(wc -l < ori.recommended.txt) markers"
echo "[INFO] ssp.snp.select.txt  : $(wc -l < ssp.snp.select.txt) markers"
echo "[INFO] adv.snp.select.txt  : $(wc -l < adv.snp.select.txt) markers"
echo

########################################
# 2. KBA version별 prefix / plink input 타입 설정
########################################

orig_mode=""
ssp_mode=""
adv_mode=""

orig_input=""
ssp_input=""
adv_input=""

orig_bim=""
ssp_bim=""
adv_bim=""

orig_extract=""
ssp_extract=""
adv_extract=""

case "${kba_version}" in
    KBAv1.1)
        orig_mode="--file"
        ssp_mode="--file"
        adv_mode="--file"

        orig_input="${ori_geno_dir}/plink/Axiom_KBAv1.1"
        ssp_input="${ssp_geno_dir}/plink/Axiom_KBAv1.1_SSP"
        adv_input="${adv_dir}/plink/Axiom_KBAv1.1_adv"

        # KBAv1.1 은 classification ID를 그대로 사용
        orig_extract="${ori_geno_dir}/classification/Recommended.ps"
        ssp_extract="ssp.snp.select.txt"
        adv_extract="adv.snp.select.txt"
        ;;
    KBAv2.0A)
        orig_mode="--bfile"
        ssp_mode="--bfile"
        adv_mode="--bfile"

        orig_input="${ori_geno_dir}/vcf/Axiom_KBAv2.0A.FINAL_genocall"
        ssp_input="${ssp_geno_dir}/vcf/Axiom_KBAv2.0A.FINAL_SSP"
        adv_input="${adv_dir}/vcf/Axiom_KBAv2.0A.FINAL_adv"

        orig_bim="${orig_input}.bim"
        ssp_bim="${ssp_input}.bim"
        adv_bim="${adv_input}.bim"

        # 실제 extract용 marker ID 생성
        orig_extract="ori.recommended.extract.txt"
        ssp_extract="ssp.snp.select.extract.txt"
        adv_extract="adv.snp.select.extract.txt"

        make_extract_from_bim "ori.recommended.txt" "$orig_bim" "$orig_extract"
        make_extract_from_bim "ssp.snp.select.txt" "$ssp_bim" "$ssp_extract"
        make_extract_from_bim "adv.snp.select.txt" "$adv_bim" "$adv_extract"
        ;;
    KBAv2.0B)
        orig_mode="--bfile"
        ssp_mode="--bfile"
        adv_mode="--bfile"

        orig_input="${ori_geno_dir}/vcf/Axiom_KBAv2.0B.FINAL_genocall"
        ssp_input="${ssp_geno_dir}/vcf/Axiom_KBAv2.0B.FINAL_SSP"
        adv_input="${adv_dir}/vcf/Axiom_KBAv2.0B.FINAL_adv"

        orig_bim="${orig_input}.bim"
        ssp_bim="${ssp_input}.bim"
        adv_bim="${adv_input}.bim"

        # 실제 extract용 marker ID 생성
        orig_extract="ori.recommended.extract.txt"
        ssp_extract="ssp.snp.select.extract.txt"
        adv_extract="adv.snp.select.extract.txt"

        make_extract_from_bim "ori.recommended.txt" "$orig_bim" "$orig_extract"
        make_extract_from_bim "ssp.snp.select.txt" "$ssp_bim" "$ssp_extract"
        make_extract_from_bim "adv.snp.select.txt" "$adv_bim" "$adv_extract"
        ;;
    *)
        echo "[ERROR] Invalid KBA version: ${kba_version}"
        echo "Available KBA_version: KBAv1.1 | KBAv2.0A | KBAv2.0B"
        exit 1
        ;;
esac

########################################
# 3. version별 output 이름
########################################

orig_out="Axiom_${kba_version}_Original_call"
ssp_out="Axiom_${kba_version}_SSP_call"
adv_out="Axiom_${kba_version}_adv_call"

echo "[INFO] Original input : ${orig_input} (${orig_mode})"
echo "[INFO] SSP input      : ${ssp_input} (${ssp_mode})"
echo "[INFO] adv input      : ${adv_input} (${adv_mode})"
echo
echo "[INFO] Original extract list : ${orig_extract} ($(wc -l < "${orig_extract}") markers)"
echo "[INFO] SSP extract list      : ${ssp_extract} ($(wc -l < "${ssp_extract}") markers)"
echo "[INFO] adv extract list      : ${adv_extract} ($(wc -l < "${adv_extract}") markers)"
echo

########################################
# 4. PLINK extract
########################################

plink ${orig_mode} "${orig_input}" \
    --extract "${orig_extract}" \
    --allow-extra-chr \
    --no-fid --no-parents --no-pheno --no-sex \
    --make-bed \
    --out "${orig_out}"

plink ${ssp_mode} "${ssp_input}" \
    --extract "${ssp_extract}" \
    --allow-extra-chr \
    --no-fid --no-parents --no-pheno --no-sex \
    --make-bed \
    --out "${ssp_out}"

plink ${adv_mode} "${adv_input}" \
    --extract "${adv_extract}" \
    --allow-extra-chr \
    --no-fid --no-parents --no-pheno --no-sex \
    --make-bed \
    --out "${adv_out}"

########################################
# 5. merge
########################################

{
    echo "${ssp_out}"
    echo "${adv_out}"
} > merge.list

plink --bfile "${orig_out}" \
    --merge-list merge.list \
    --allow-extra-chr \
    --make-bed \
    --out "${merge_prefix}"

echo
echo "[DONE] Output: ${merge_prefix}.bed / .bim / .fam"