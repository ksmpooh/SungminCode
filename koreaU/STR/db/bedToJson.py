### TGRT bed file to EH json


import json,os,glob

inBED = "pathogenic_repeats.hg38.bed"
outjson = inBED.replace(".bed","bedtojson.json")


## eh.pathogenic.txt
#ID=AFF2;MOTIFS=GCC;STRUC=(GCC)n
#ID=AR;MOTIFS=GCA;STRUC=(GCA)n
#ID=ATN1;MOTIFS=CAG;STRUC=(CAG)n

refin = open("/BDATA/smkim/STR/EH/resources/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed","r")
ref = [s.split(";")[0].split("=")[1] for s in refin]

target_gene = []

def bed_to_json(bed_file, json_file):
    json_data = []

    with open(bed_file, 'r') as bed:
        for line in bed:
            fields = line.strip().split()
            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            id_field = fields[3].split(';')
            
            locus_id = id_field[0].split('=')[1]  # Extract ID
            if locus_id in ref:
                continue
            motifs = id_field[1].split('=')[1]    # Extract MOTIFS
            repeat_units = motifs.split(',')      # Support multiple motifs
            reference_region = f"{chrom}:{start}-{end}"

            if locus_id not in target_gene:
                target_gene.append(locus_id)
            else:
                locus_id = locus_id + "_" + repeat_units
            
            for repeat_unit in repeat_units:
                entry = {
                    "LocusId": locus_id,
                    "LocusStructure": repeat_unit,
                    "ReferenceRegion": reference_region,
                    "VariantType": "Repeat",
                    #"OffTargetRegions": []
                }
                json_data.append(entry)

    with open(json_file, 'w') as outfile:
        json.dump(json_data, outfile, indent=4)

bed_to_json(inBED, outjson)


for i in target_gene:
    print(i)

'''
{
        "LocusId": "AFF2",
        "LocusStructure": "(GCC)*",
        "ReferenceRegion": "chrX:148500631-148500691",
        "VariantType": "Repeat"
    },
    {
        "LocusId": "CNBP",
        "LocusStructure": "(CAGG)*(CAGA)*(CA)*",
        "ReferenceRegion": [
            "chr3:129172576-129172656",
            "chr3:129172656-129172696",
            "chr3:129172696-129172732"
        ],
        "VariantId": [
            "CNBP",
            "CNBP_CAGA",
            "CNBP_CA"
        ],
        "VariantType": [
            "Repeat",
            "Repeat",
            "Repeat"
        ]
    },
'''
'''
chr3    138946020       138946063       ID=FOXL2;MOTIFS=NGC;STRUC=(NGC)n
chr3    183712187       183712223       ID=YEATS2;MOTIFS=TTTTA,TTTCA;STRUC=(TTTTA)n(TTTCA)n
chr4    3074876 3074966 ID=HTT;MOTIFS=CAG,CCG;STRUC=(CAG)nCAACAG(CCG)n
chr4    39348424        39348479        ID=RFC1;MOTIFS=AAAAG,AAAGG,AAGGG,AAGAG,AGAGG,AACGG,GGGAC,AAAGGG;STRUC=<RFC1>
chr4    41745972        41746032        ID=PHOX2B;MOTIFS=GCN;STRUC=(GCN)n
chr4    159342526       159342617       ID=RAPGEF2;MOTIFS=TTTTA,TTTCA;STRUC=(TTTTA)n(TTTCA)n(TTTTA)n
chr5    10356346        10356412        ID=MARCHF6;MOTIFS=TTTTA,TTTCA;STRUC=(TTTTA)n(TTTCA)n

chr3	63912684	63912726	ID=ATXN7;MOTIFS=GCA,GCC;STRUC=(GCA)n(GCC)n	ATXN7	(GCA)*(GCC)+	chr3:63912684-63912714_chr3:63912714-63912726
chr3	129172576	129172732	ID=CNBP;MOTIFS=CAGG,CAGA,CA;STRUC=(CAGG)n(CAGA)n(CA)n	CNBP	(CAGG)*(CAGA)*(CA)*	chr3:129172576-129172656_chr3:129172656-129172696_chr3:129172696-129172732
chr9	69037261	69037304	ID=FXN;MOTIFS=A,GAA;STRUC=(A)n(GAA)n	FXN	(A)*(GAA)*	chr9:69037261-69037286_chr9:69037286-69037304
chr4	3074876	3074966	ID=HTT;MOTIFS=CAG,CCG;STRUC=(CAG)nCAACAG(CCG)n	HTT	(CAG)*CAACAG(CCG)*	chr4:3074876-3074933_chr4:3074939-3074966
chr20	2652733	2652775	ID=NOP56;MOTIFS=GGCCTG,CGCCTG;STRUC=(GGCCTG)n(CGCCTG)n	NOP56	(GGCCTG)*(CGCCTG)*	chr20:2652733-2652757_chr20:2652757-2652775
'''

###
#(base) ➜  db head eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed
#chrX	148500631	148500691	ID=AFF2;MOTIFS=GCC;STRUC=(GCC)n
#chrX	67545316	67545385	ID=AR;MOTIFS=GCA;STRUC=(GCA)n
#chr12	6936716	6936773	ID=ATN1;MOTIFS=CAG;STRUC=(CAG)n

