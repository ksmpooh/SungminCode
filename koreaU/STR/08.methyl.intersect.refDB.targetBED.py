import os,glob

file_a = 'Fianl.ref.sort.STR.prop.bed'  # First BED file
file_b = 'NIH23F1013274.pbmm2_hg38_withunmapped.pb_CpG.combined.bed'  # Second BED file
output_file = 'test.bed'


# Define the function to parse BED files
def read_bed_file(filepath):
    with open(filepath, 'r') as f:
        return [line.strip().split('\t') for line in f.readlines()]

# Function to find intersections
def intersect_bed_files(bed_a, bed_b, output_path):
    for chrom_a, start_a, end_a, *rest_a in bed_a:
        start_a, end_a = int(start_a), int(end_a)
        for chrom_b, start_b, end_b, *rest_b in bed_b:
            start_b, end_b = int(start_b), int(end_b)
                # Check if ranges overlap and are on the same chromosome
            if chrom_a == chrom_b and start_b < end_a and end_b > start_a:                   
                out = '\t'.join([chrom_b, str(max(start_a, start_b)), str(min(end_a, end_b))] + rest_b) + '\n'
                os.system("%s >> %s"%(out,output_path) )
                #output.write('\t'.join([chrom_b, str(max(start_a, start_b)), str(min(end_a, end_b))] + rest_b) + '\n')

# Paths to input files

# Read the BED files
bed_a = read_bed_file(file_a)
bed_b = read_bed_file(file_b)

# Perform the intersection and write the output
intersect_bed_files(bed_a, bed_b, output_file)