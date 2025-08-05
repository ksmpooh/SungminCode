#!/bin/bash
#ls *polished | xargs -I {} -P 80 bash -c "bash liftoff.after.sh {} CHM13 or GRCh38"


input=$1       # GFF3 input file
build=$2       # CHM13 또는 GRCh38

mkdir -p liftoff_processing/transcript/
mkdir -p liftoff_processing/gene/

# --- Transcript 블럭 ---
if [ "$build" == "CHM13" ]; then
  awk '$3=="transcript" {
    attr = $NF
    source_gene = "NA"
    transcript_biotype = "NA"
    source_transcript = "NA"

    n = split(attr, fields, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", fields[i])
        if (fields[i] ~ /^source_gene=/) {
            split(fields[i], a, "="); source_gene = a[2]
        } else if (fields[i] ~ /^gene_biotype=/) {
            split(fields[i], a, "="); transcript_biotype = a[2]
        } else if (fields[i] ~ /^source_transcript=/) {
            split(fields[i], a, "="); source_transcript = a[2]
        }
        
    }
    print $1, source_gene, transcript_biotype, source_transcript
  }' "$input" > liftoff_processing/transcript/$(basename "$input").afterliftoff_transcript

else  # GRCh38
  awk '$3=="transcript" {
    attr = $NF
    parent = "NA"
    transcript_type = "NA"
    transcript_id = "NA"
    level = "NA"

    n = split(attr, fields, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", fields[i])
        if (fields[i] ~ /^Parent=/) {
            split(fields[i], a, "="); parent = a[2]
        } else if (fields[i] ~ /^transcript_type=/) {
            split(fields[i], a, "="); transcript_type = a[2]
        } else if (fields[i] ~ /^transcript_id=/) {
            split(fields[i], a, "="); transcript_id = a[2]
        } else if (fields[i] ~ /^level=/) {
            split(fields[i], a, "="); level = a[2]
        }
    }
    print $1, parent, transcript_type, transcript_id,level
  }' "$input" > liftoff_processing/transcript/$(basename "$input").afterliftoff_transcript
fi

# --- Gene 블럭 ---
if [ "$build" == "CHM13" ]; then
  awk '$3=="gene" {
    attr = $NF
    gene_name = "NA"
    gene_biotype = "NA"
    extra_copy_number = "NA"
    coverage = "NA"

    n = split(attr, fields, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", fields[i])
        if (fields[i] ~ /^gene_name=/) {
            split(fields[i], a, "="); gene_name = a[2]
        } else if (fields[i] ~ /^gene_biotype=/) {
            split(fields[i], a, "="); gene_biotype = a[2]
        } else if (fields[i] ~ /^extra_copy_number=/) {
            split(fields[i], a, "="); extra_copy_number = a[2]
        } else if (fields[i] ~ /^coverage=/) {
            split(fields[i], a, "="); coverage = a[2]
        }
    }
    print $1, gene_name, gene_biotype, extra_copy_number, coverage
  }' "$input" > liftoff_processing/gene/$(basename "$input").afterliftoff_gene

else  # GRCh38
  awk '$3=="gene" {
    attr = $NF
    gene_name = "NA"
    gene_biotype = "NA"
    extra_copy_number = "NA"
    coverage = "NA"
    level = "NA"
    valid_ORFs = "NA"

    n = split(attr, fields, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", fields[i])
        if (fields[i] ~ /^gene_name=/) {
            split(fields[i], a, "="); gene_name = a[2]
        } else if (fields[i] ~ /^gene_type=/) {
            split(fields[i], a, "="); gene_biotype = a[2]
        } else if (fields[i] ~ /^extra_copy_number=/) {
            split(fields[i], a, "="); extra_copy_number = a[2]
        } else if (fields[i] ~ /^coverage=/) {
            split(fields[i], a, "="); coverage = a[2]
        }else if (fields[i] ~ /^level=/) {
            split(fields[i], a, "="); level = a[2]
        }else if (fields[i] ~ /^valid_ORFs=/) {
            split(fields[i], a, "="); valid_ORFs = a[2]
        }
    }
    print $1, gene_name, gene_biotype, extra_copy_number, coverage, level, valid_ORFs
  }' "$input" > liftoff_processing/gene/$(basename "$input").afterliftoff_gene
fi