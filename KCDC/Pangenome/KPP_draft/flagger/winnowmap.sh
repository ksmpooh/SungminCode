#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# End-to-end: meryl -> repetitive kmers -> winnowmap -> sorted BAM
#
# Inputs:
#   1) assembly.fa
#   2) reads.fastq.gz
#   3) output_prefix
#
# Example:
#   bash run_winnowmap_with_meryl.sh asm.fa reads.fastq.gz KPPD132
###############################################################################

ASSEMBLY_FA="$1"
READS_FASTQ="$2"
OUT_PREFIX="$3"

# ------------------ tools ------------------
MERYL="meryl"
WINNOWMAP="winnowmap"
SAMTOOLS="samtools"

# ------------------ parameters ------------------
kmerSize=15
merylThreads=32
merylMemory=128000        # MB (128 GB)
distinctCutoff=0.9998

preset="map-pb"
options="--eqx --cs -Y -L -y -I8g -p0.5"
threads=100

# ------------------ outputs ------------------
MERYL_DB="${OUT_PREFIX}.k${kmerSize}.merylDB"
REPETITIVE_TXT="${OUT_PREFIX}.repetitive_k${kmerSize}.txt"
SORTED_BAM="${OUT_PREFIX}.sorted.bam"

###############################################################################
# 1) meryl count (multithread + memory)
###############################################################################
echo "[STEP 1] meryl count (k=${kmerSize})"

"${MERYL}" count \
  "k=${kmerSize}" \
  threads=${merylThreads} \
  memory=${merylMemory} \
  output "${MERYL_DB}" \
  "${ASSEMBLY_FA}"

###############################################################################
# 2) extract repetitive kmers
###############################################################################
echo "[STEP 2] meryl print (distinct > ${distinctCutoff})"

"${MERYL}" print greater-than distinct=${distinctCutoff} \
  "${MERYL_DB}" \
  > "${REPETITIVE_TXT}"

###############################################################################
# 3) winnowmap -> BAM -> sort
###############################################################################
echo "[STEP 3] winnowmap + samtools sort"

"${WINNOWMAP}" \
  -W "${REPETITIVE_TXT}" \
  -a \
  -x "${preset}" \
  ${options} \
  -t${threads} \
  "${ASSEMBLY_FA}" \
  "${READS_FASTQ}" \
| "${SAMTOOLS}" view -h -b - \
| "${SAMTOOLS}" sort -@${threads} -o "${SORTED_BAM}" -

###############################################################################
# 4) (Optional) index
###############################################################################
"${SAMTOOLS}" index -@${threads} "${SORTED_BAM}"

echo "[DONE]"
echo "  BAM:  ${SORTED_BAM}"
echo "  BAI:  ${SORTED_BAM}.bai"
echo "  kmers:${REPETITIVE_TXT}"