#pbmm2 align reference.fasta input.bam output.bam --sort --preset CCS -L 0.1 -c 0

#!/bin/bash
# ---------------------------------------------------
# DeepVariant pipeline (fixed output naming)
# Author: smkim style
#
# Usage:
# bash run_deepvariant_pipeline.sh <ref.fasta> <input.bam> <output_prefix> <sample_id> <theme>
#
# Example:
# bash run_deepvariant_pipeline.sh \
#   /data/ref/HLA.target.fasta \
#   /data/reads/KPPD019.sorted.bam \
#   /data/out/KPPD019.short_read.Deepvariant \
#   KPPD019 \
#   short
#
# Theme:
#   short -> model_type=WGS
#   long  -> model_type=PACBIO
#
# Output files:
#   <output_prefix>.g.vcf.gz
#   <output_prefix>.vcf.gz
# ---------------------------------------------------

set -euo pipefail

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <ref.fasta> <input.bam> <output_prefix> <sample_id> <theme: short|long>"
  exit 1
fi

REF_FA=$1
INPUT_BAM=$2
OUT_PREFIX=$3   # e.g. /outDir/sampleID.short_read.Deepvariant
SAMPLE_ID=$4
THEME=$5         # short or long

NUM_SHARDS=${NUM_SHARDS:-16}
DOCKER_IMAGE="google/deepvariant:1.9.0"

# -------------------------------
# path resolution
# -------------------------------
REF_DIR=$(dirname "${REF_FA}")
INPUT_DIR=$(dirname "${INPUT_BAM}")
OUT_DIR=$(dirname "${OUT_PREFIX}")
BASE_NAME=$(basename "${OUT_PREFIX}")

OUT_VCF="${OUT_PREFIX}.vcf.gz"
OUT_GVCF="${OUT_PREFIX}.g.vcf.gz"
LOGFILE="${OUT_PREFIX}.log"

echo "---------------------------------------------"
echo "[DEEPVARIANT START]"
echo "🧬 Reference : ${REF_FA}"
echo "📦 BAM        : ${INPUT_BAM}"
echo "📤 Out prefix : ${OUT_PREFIX}"
echo "🧩 Sample ID  : ${SAMPLE_ID}"
echo "🎨 Theme      : ${THEME}"
echo "🧵 Threads    : ${NUM_SHARDS}"
echo "🐳 Docker     : ${DOCKER_IMAGE}"
echo "---------------------------------------------"

# -------------------------------
# Index check
# -------------------------------
if [ ! -f "${REF_FA}.fai" ]; then
  echo "[INFO] Creating FASTA index..."
  samtools faidx "${REF_FA}"
fi

if [ ! -f "${INPUT_BAM}.bai" ] && [ ! -f "${INPUT_BAM%.bam}.bai" ]; then
  echo "[INFO] Creating BAM index..."
  samtools index "${INPUT_BAM}"
fi

# -------------------------------
# Model type
# -------------------------------
case "${THEME,,}" in
  short)
    MODEL_TYPE="WGS"
    ;;
  long)
    MODEL_TYPE="PACBIO"
    ;;
  *)
    echo "Theme must be 'short' or 'long'."
    exit 1
    ;;
esac

echo "[INFO] Using model_type=${MODEL_TYPE}"

# -------------------------------
# Run DeepVariant (Docker)
# -------------------------------
CONTAINER_REF="/ref/$(basename "${REF_FA}")"
CONTAINER_READS="/input/$(basename "${INPUT_BAM}")"
CONTAINER_OUT_VCF="/output/${BASE_NAME}.vcf.gz"
CONTAINER_OUT_GVCF="/output/${BASE_NAME}.g.vcf.gz"

docker run --rm \
  -v "${INPUT_DIR}":/input:ro \
  -v "${OUT_DIR}":/output \
  -v "${REF_DIR}":/ref:ro \
  "${DOCKER_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type "${MODEL_TYPE}" \
    --ref "${CONTAINER_REF}" \
    --reads "${CONTAINER_READS}" \
    --output_vcf "${CONTAINER_OUT_VCF}" \
    --output_gvcf "${CONTAINER_OUT_GVCF}" \
    --num_shards "${NUM_SHARDS}" \
  2>&1 | tee "${LOGFILE}"

# -------------------------------
# Index VCFs
# -------------------------------
if [ -f "${OUT_VCF}" ]; then
  echo "[INFO] Indexing ${OUT_VCF}..."
  tabix -p vcf "${OUT_VCF}"
fi

if [ -f "${OUT_GVCF}" ]; then
  echo "[INFO] Indexing ${OUT_GVCF}..."
  tabix -p vcf "${OUT_GVCF}"
fi

echo "---------------------------------------------"
echo "[✅ DONE] DeepVariant completed"
echo "Output:"
echo "  ${OUT_VCF}"
echo "  ${OUT_GVCF}"
echo "Log:"
echo "  ${LOGFILE}"
echo "---------------------------------------------"