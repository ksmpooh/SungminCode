#!/usr/bin/env bash
# -----------------------------------------
# One-pass alignment & GATK duplicate marking (minimal temp)
# Author: you
# Usage:
#   bash run_gatk_align_dedup.sh <ref.fa> <reads_R1.fq[.gz]> <reads_R2.fq[.gz]> <sample_id> <out_dir> [threads]
#
# Output:
#   <out_dir>/<sample_id>.dedup.bam
#   <out_dir>/<sample_id>.dedup.bam.bai
#   <out_dir>/<sample_id>.dup_metrics.txt
# -----------------------------------------
set -euo pipefail

REF=$1          # reference FASTA (e.g., chm13v2.0.fa)
R1=$2           # FASTQ R1 (gz/fastq ok)
R2=$3           # FASTQ R2 (gz/fastq ok)
SAMPLE=$4       # sample id (e.g., KPPD019)
OUTDIR=$5       # output directory
THREADS=${6:-32}

mkdir -p "${OUTDIR}"
LOG="${OUTDIR}/${SAMPLE}.log"
exec > >(tee -a "${LOG}") 2>&1

echo "---------------------------------------------"
echo "[GATK PIPELINE START]"
echo "🧬 REF     : ${REF}"
echo "📄 R1      : ${R1}"
echo "📄 R2      : ${R2}"
echo "📛 SAMPLE  : ${SAMPLE}"
echo "📂 OUTDIR  : ${OUTDIR}"
echo "🧵 THREADS : ${THREADS}"
echo "---------------------------------------------"

# --------- Preflight: reference indexes ----------
# 1) bwa index
if ! [ -f "${REF}.bwt" ] && ! [ -f "${REF}.0123" ]; then
  echo "[INFO] bwa index not found. Building..."
  bwa index "${REF}"
fi

# 2) samtools faidx
if ! [ -f "${REF}.fai" ]; then
  echo "[INFO] faidx not found. Building..."
  samtools faidx "${REF}"
fi

# 3) GATK/Picard sequence dictionary
DICT="${REF%.*}.dict"
if ! [ -f "${DICT}" ]; then
  echo "[INFO] sequence dictionary not found. Building..."
  gatk CreateSequenceDictionary -R "${REF}" -O "${DICT}"
fi

# --------- Temp space (auto-clean) ----------
TMPBASE="${OUTDIR}/.tmp_${SAMPLE}_$(date +%s)"
mkdir -p "${TMPBASE}"
cleanup() {
  echo "[CLEANUP] removing temp directory ${TMPBASE}"
  rm -rf "${TMPBASE}" || true
}
trap cleanup EXIT

# temp sorted BAM (intermediate) — will be deleted at the end by trap
SORTED_BAM="${TMPBASE}/${SAMPLE}.sorted.bam"

# Final outputs
DEDUP_BAM="${OUTDIR}/${SAMPLE}.dedup.bam"
DEDUP_METRICS="${OUTDIR}/${SAMPLE}.dup_metrics.txt"

# --------- Read Group ----------
# (필요 시 LB/PU 수정)
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPU:${SAMPLE}.1"

# --------- Step 1: Align & sort (streaming) ----------
# bwa mem -> samtools sort (no on-disk SAM)
echo "[STEP 1] bwa mem → samtools sort"
bwa mem -t "${THREADS}" -R "${RG}" "${REF}" "${R1}" "${R2}" \
  | samtools sort -@ "${THREADS}" -m 3G -o "${SORTED_BAM}" -

# --------- Step 2: MarkDuplicates (GATK) ----------
# GATK MarkDuplicatesSpark expects a file; we use the sorted BAM above and remove it afterwards.
echo "[STEP 2] GATK MarkDuplicatesSpark"
gatk --java-options "-Xmx16g -Djava.io.tmpdir=${TMPBASE}" MarkDuplicatesSpark \
  -I "${SORTED_BAM}" \
  -O "${DEDUP_BAM}" \
  -M "${DEDUP_METRICS}" \
  --create-output-bam-index true \
  --tmp-dir "${TMPBASE}" \
  --spark-master "local[${THREADS}]"

# Optional: remove the intermediate immediately (trap will also remove it)
rm -f "${SORTED_BAM}" || true

echo "---------------------------------------------"
echo "[✅ DONE] Outputs:"
echo "  BAM : ${DEDUP_BAM}"
echo "  BAI : ${DEDUP_BAM}.bai"
echo "  MET : ${DEDUP_METRICS}"
echo "  LOG : ${LOG}"
echo "---------------------------------------------"