#!/bin/bash
# -------------------------------
# Flagger annotation (post-analysis) pipeline
# Author: smkim
# Description:
#   Reuse existing flagger mapping results (sorted.bam) and run CHM13-based
#   annotation liftover, bias detection (bam2cov), and HMM-Flagger.
#
# Usage:
#   bash run_flagger_annotation_pipeline.sh <assembly.fa> <ref_dir> <flagger_output_dir> <sample_id>
#
# Example:
#   bash run_flagger_annotation_pipeline.sh \
#     /CDATA/pangenome/flagger_test/assembly/KPPD019.assembly.merge.fa \
#     /CDATA/pangenome/flagger_test/ref \
#     /CDATA/pangenome/flagger_test/flagger_output \
#     KPPD019
# -------------------------------

set -euo pipefail

# ---------- Input arguments ----------
ASSEMBLY_FA=$1      # Assembly fasta
REF_DIR=$2          # Directory containing alpha TSV + CHM13 annotations
OUT_BASE=$3         # Directory containing previous flagger outputs
SAMPLE_ID=$4        # e.g. KPPD019

# ---------- Derived paths ----------
ASSEMBLY_DIR=$(dirname "${ASSEMBLY_FA}")
FLAGGER_OUT="${OUT_BASE}/${SAMPLE_ID}_flagger_out"       # Existing flagger result folder
ANNOT_DIR="${OUT_BASE}/${SAMPLE_ID}_annotation"           # New annotation analysis folder
LIFTED_DIR="${ANNOT_DIR}/lifted"
mkdir -p "${LIFTED_DIR}"

BAM="${FLAGGER_OUT}/${SAMPLE_ID}.sorted.bam"
BED="${FLAGGER_OUT}/${SAMPLE_ID}.bed"
COV="${ANNOT_DIR}/${SAMPLE_ID}.cov.gz"
JSON="${ANNOT_DIR}/${SAMPLE_ID}.json"
PAF="${ANNOT_DIR}/${SAMPLE_ID}.ref2asm_corrected.paf"
ALPHA_TSV="${REF_DIR}/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv"

echo "---------------------------------------------"
echo "[FLAGGER ANNOTATION PIPELINE START]"
echo "📘 Assembly fasta  : ${ASSEMBLY_FA}"
echo "📄 BAM file        : ${BAM}"
echo "📄 BED file        : ${BED}"
echo "📂 Reference dir   : ${REF_DIR}"
echo "📤 Output dir      : ${ANNOT_DIR}"
echo "📛 Sample ID       : ${SAMPLE_ID}"
echo "---------------------------------------------"
sleep 1

# ---------- Step 1. Make PAF (CHM13 → assembly) ----------
echo "[INFO] Generating CHM13→assembly PAF..."
minimap2 -x asm20 -c -t 32 "${ASSEMBLY_FA}" "${REF_DIR}/chm13v2.0_maskedY_rCRS.fa" > "${PAF}"

# ---------- Step 2. LiftOver CHM13 annotations ----------
echo "[INFO] Lifting CHM13 annotations..."
for bed in \
  "${REF_DIR}/potential_biases/chm13v2.0_bsat.bed" \
  "${REF_DIR}/potential_biases/chm13v2.0_hsat1A.bed" \
  "${REF_DIR}/potential_biases/chm13v2.0_hsat1B.bed" \
  "${REF_DIR}/potential_biases/chm13v2.0_hsat2.bed" \
  "${REF_DIR}/potential_biases/chm13v2.0_hsat3.bed" \
  "${REF_DIR}/potential_biases/chm13v2.0_hor.bed" \
  "${REF_DIR}/stratifications/censat/chm13v2.0_censat.bed" \
  "${REF_DIR}/stratifications/sd/chm13v2.0_SD.all.bed" \
  "${REF_DIR}/stratifications/sex/chm13v2.0_sex.bed" \
  "${REF_DIR}/stratifications/repeat_masker/chm13v2.0_RM_4.1.2p1_le6_STR.bed" \
  "${REF_DIR}/stratifications/repeat_masker/chm13v2.0_RM_4.1.2p1_ge7_VNTR.bed"
do
  base=$(basename "$bed" .bed)
  echo "   lifting $base ..."
  paftools.js liftover "${PAF}" "$bed" > "${LIFTED_DIR}/${base}.lifted.bed"
done

# ---------- Step 3. Build JSON for bam2cov ----------
echo "[INFO] Creating annotation JSON..."
cat > "${JSON}" <<EOF
{
  "whole_genome": "/wdir/${SAMPLE_ID}_flagger_out/${SAMPLE_ID}.bed",
  "bsat": "/ref_lifted/chm13v2.0_bsat.lifted.bed",
  "hsat1A": "/ref_lifted/chm13v2.0_hsat1A.lifted.bed",
  "hsat1B": "/ref_lifted/chm13v2.0_hsat1B.lifted.bed",
  "hsat2": "/ref_lifted/chm13v2.0_hsat2.lifted.bed",
  "hsat3": "/ref_lifted/chm13v2.0_hsat3.lifted.bed",
  "hor": "/ref_lifted/chm13v2.0_hor.lifted.bed",
  "censat": "/ref_lifted/chm13v2.0_censat.lifted.bed",
  "sd": "/ref_lifted/chm13v2.0_SD.all.lifted.bed",
  "sex": "/ref_lifted/chm13v2.0_sex.lifted.bed",
  "VNTR": "/ref_lifted/chm13v2.0_RM_4.1.2p1_ge7_VNTR.lifted.bed",
  "STR": "/ref_lifted/chm13v2.0_RM_4.1.2p1_le6_STR.lifted.bed"
}
EOF

# ---------- Step 4. Run bam2cov (bias detection) ----------
echo "[INFO] Running bam2cov with bias detection..."
docker run -it --rm \
  -v "${OUT_BASE}":/wdir \
  -v "${LIFTED_DIR}":/ref_lifted \
  -v "${ANNOT_DIR}":/outdir \
  mobinasri/flagger:v1.1.0 \
  bam2cov \
    --bam /wdir/${SAMPLE_ID}_flagger_out/${SAMPLE_ID}.sorted.bam \
    --output /outdir/${SAMPLE_ID}.cov.gz \
    --threads 48 \
    --annotationJson /outdir/${SAMPLE_ID}.json \
    --baselineAnnotation whole_genome \
    --runBiasDetection

# ---------- Step 5. Run HMM Flagger ----------
echo "[INFO] Running HMM flagger..."
mkdir -p "${ANNOT_DIR}/final_results"

docker run -it --rm \
  -v "${ANNOT_DIR}":/work \
  -v "${REF_DIR}":/ref \
  mobinasri/flagger:v1.1.0 \
  hmm_flagger \
    --input /work/${SAMPLE_ID}.cov.gz \
    --outputDir /work/final_results \
    --alphaTsv /ref/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv \
    --labelNames Err,Dup,Hap,Col \
    --modelType trunc_exp_gaussian \
    --threads 48

echo "---------------------------------------------"
echo "[✅ DONE] Results saved in: ${ANNOT_DIR}/final_results"
echo "---------------------------------------------"


#    --minHighMapqRatio 0.5 \
#    --maxHighMapqRatio 0.25 \