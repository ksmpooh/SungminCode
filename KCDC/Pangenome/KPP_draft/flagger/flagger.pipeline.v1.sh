#!/bin/bash
# -------------------------------
# Flagger HMM pipeline (fixed output name)
# Author: smkim
# Usage:
# bash run_flagger_pipeline.sh <assembly.fa> <reads.fastq.gz> <ref_dir> <flagger_output_dir> <sample_id>
# Example:
# bash run_flagger_pipeline.sh \
#   /CDATA/pangenome/flagger_test/assembly/KPPD019.assembly.merge.fa \
#   /CDATA/pangenome/flagger_test/raw_fastq/NIH23F1013274.Revio_hifi_reads.rawfastq.merge.fastq.gz \
#   /CDATA/pangenome/flagger_test/ref \
#   /CDATA/pangenome/flagger_test/flagger_output \
#   KPPD019
# -------------------------------

set -euo pipefail

# ---------- Input args ----------
ASSEMBLY_FA=$1       # absolute path to assembly fasta
READS_FASTQ=$2       # absolute path to reads fastq
REF_DIR=$3           # directory containing alpha TSV
OUT_BASE=$4          # base output directory
SAMPLE_ID=$5         # sample name (used for output folder)

# ---------- Paths ----------
ASSEMBLY_DIR=$(dirname "${ASSEMBLY_FA}")
READS_DIR=$(dirname "${READS_FASTQ}")
ALPHA_TSV="${REF_DIR}/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv"
OUTDIR="${OUT_BASE}/${SAMPLE_ID}_flagger_out"

mkdir -p "${OUTDIR}"

echo "---------------------------------------------"
echo "[FLAGGER PIPELINE START]"
echo "📘 Assembly fasta : ${ASSEMBLY_FA}"
echo "📂 Assembly dir   : ${ASSEMBLY_DIR}"
echo "📄 Reads file     : ${READS_FASTQ}"
echo "📂 Reads dir      : ${READS_DIR}"
echo "📄 Alpha TSV      : ${ALPHA_TSV}"
echo "📤 Output dir     : ${OUTDIR}"
echo "📛 Sample ID      : ${SAMPLE_ID}"
echo "---------------------------------------------"
sleep 1

# ---------- Index & BED ----------
FAI="${ASSEMBLY_FA}.fai"
BED="${OUTDIR}/${SAMPLE_ID}.bed"

if [ ! -f "${FAI}" ]; then
  echo "[INFO] Generating fasta index..."
  samtools faidx "${ASSEMBLY_FA}"
fi

awk '{print $1 "\t0\t" $2}' "${FAI}" > "${BED}"

# ---------- Mapping ----------
echo "[INFO] Mapping HiFi reads to assembly..."
BAM="${OUTDIR}/${SAMPLE_ID}.sorted.bam"

minimap2 --cs -L -Y -t 64 -ax map-hifi "${ASSEMBLY_FA}" "${READS_FASTQ}" \
  | samtools sort -@ 16 -o "${BAM}"
samtools index "${BAM}"

# ---------- Coverage ----------
echo "[INFO] Generating coverage..."
JSON="${OUTDIR}/annotations_path.json"
COV="${OUTDIR}/${SAMPLE_ID}.cov.gz"

cat > "${JSON}" <<EOF
{
  "whole_genome": "/work/${SAMPLE_ID}.bed"
}
EOF

docker run -it --rm \
  -v "${OUTDIR}":/work \
  mobinasri/flagger:v1.1.0 \
  bam2cov \
    --bam /work/${SAMPLE_ID}.sorted.bam \
    --output /work/${SAMPLE_ID}.cov.gz \
    --threads 48 \
    --annotationJson /work/annotations_path.json \
    --baselineAnnotation whole_genome

# ---------- HMM Flagger ----------
echo "[INFO] Running HMM flagger..."


mkdir -p "${OUTDIR}"/final_results

docker run -it --rm \
  -v "${OUTDIR}":/work \
  -v "${REF_DIR}":/ref \
  mobinasri/flagger:v1.1.0 \
  hmm_flagger \
    --input /work/${SAMPLE_ID}.cov.gz \
    --outputDir /work/final_results \
    --alphaTsv /ref/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv \
    --labelNames Err,Dup,Hap,Col \
    --threads 48

echo "---------------------------------------------"
echo "[✅ DONE] Results saved in: ${OUTDIR}/final_results"
echo "---------------------------------------------"