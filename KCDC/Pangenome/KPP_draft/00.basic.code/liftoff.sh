#conda activate liftoff_env
#ls *fa |sed 's/.fa//g' | xargs -I{} -P 6 bash -c "bash liftoff.sh {} CHM13 outDir"
# memery 50G per job

# Parse arguments
build=$2  # CHM13 or GRCh38
outDir=$3  # Output directory

# Set the reference genome and GFF file based on the build
if [ "$build" == "GRCh38" ]; then
    REFERENCE=/ADATA/pangenome/db/gencode/GRCh38.primary_assembly.genome.fa
    theme="gencode_GRCh38"
    GFF_gencode=/ADATA/pangenome/db/gencode/gencode.v47.primary_assembly.annotation.gff3
    echo "Using GRCh38 reference genome and GFF file"
elif [ "$build" == "CHM13" ]; then
    REFERENCE=/ADATA/pangenome/db/chm13/chm13v2.0_maskedY_rCRS.fa
    theme="gencode_CHM13"
    GFF_gencode=/ADATA/pangenome/db/chm13/db_test/chm13.draft_v2.0.gene_annotation.gff3_db
    echo "Using T2TCHM13 reference genome and GFF file"
else
    echo "Invalid theme specified. Please use either 'CHM13' or 'GRCh38'."
    exit 1
fi

# Create the output directory for the theme
mkdir -p $outDir/$theme

# Set file names and paths
TARGET=$1.fa
OUTPUT=$outDir/$theme/$1.$theme.liftoff.gff3
UNMAP=$outDir/$theme/$1.$theme.unMap.liftoff.txt
THREADS=16

# Run Liftoff
liftoff -sc 0.90 -copies -g $GFF_gencode -o $OUTPUT -u $UNMAP -p $THREADS -dir $outDir/$theme/$1 -polish $TARGET $REFERENCE

# Clean up temporary files (optional)
rm -rf $outDir/$theme/$1