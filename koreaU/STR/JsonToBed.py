### STR DB format change: JSON (EH) to BED (TRGT) (2024-05-31 by smkim)
### python jsonToBed.py [input.JSON]


import os, glob, sys
import json
import re


#wDir= "/BDATA/smkim/STR/EH/resources/"
#datain = wDir + "eh.v5_w_gangstr.v13.polymorphic.json"
datain = sys.argv[1]
with open(datain, 'r') as f:
    json_data = json.load(f)

def jsontobed_byline_forcheck(i):
    LocusId = i["LocusId"]
    LocusStructure = i["LocusStructure"]
    ReferenceRegion = i["ReferenceRegion"]
    print(LocusId,LocusStructure,ReferenceRegion)
    if type(ReferenceRegion) == list:
        ReferenceRegion = '_'.join(ReferenceRegion)
    #chrom, positions = ReferenceRegion.split(':')
    #start, end = positions.split('-')
    tmp = [LocusId,LocusStructure,ReferenceRegion]
#    print("tmp")
    print(tmp)
    return tmp

def data_check():
    out = open("Common_trgt_repeat_add_eh.v5_w_gangstr.v13.polymorphic.txt","w")
    tmp = []
    bed = open("/BDATA/smkim/TOOLs/trgt/repeats/pathogenic_repeats.hg38.bed","r")
    bed = bed.readlines()
    idlist = [s.split()[3].split(";")[0].split('=')[1] for s in bed]
    bed_dic = {}
    for i in bed:
        bed_dic[i.split()[3].split(";")[0].split('=')[1]]=i.strip()
    #print(bed_dic)
    for i in json_data:
        if i['LocusId'] in idlist:
            out.write('%s\t%s\n'%(bed_dic[i['LocusId']],'\t'.join(jsontobed_byline_forcheck(i))))
            tmp.append(i)
    out.close()
    return tmp


## EH vs TRGT 공통된 STR 를 추출하고 비교하여 패턴을 확인 후 변경
gene_list = data_check()

''' 양쪽 DB 모두에 있는 STR ID
AFF2 AR ATN1 ATXN1 ATXN10 ATXN2 ATXN3 ATXN7 C9ORF72 CACNA1A CBL CNBP CSTB DIP2B DMPK FMR1 FXN GIPC1 GLS HTT JPH3 NIPA1 NOP56 PABPN1 PHOX2B PPP2R2B RFC1 TBP TCF4
## check 
ATXN7 C9ORF72 CNBP FMR1 FXN HTT NOP56
grep -E ATXN7|C9ORF72|CNBP|FMR1|FXN|HTT|NOP56
'''
#pathogenic_repeats.hg38.bed:chrX	148500631	148500691	ID=AFF2;MOTIFS=GCC;STRUC=(GCC)n
#pathogenic_repeats.hg38.bed:chr3	129172576	129172732	ID=CNBP;MOTIFS=CAGG,CAGA,CA;STRUC=(CAGG)n(CAGA)n(CA)n
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
## extract motif from STRUC
def extract_motif(text):
    # Regular expression to find all substrings inside parentheses
    pattern = r'\(([^)]+)\)'
    # Find all matching substrings
    matches = re.findall(pattern, text)
    # Join matches with a comma separator
    result = ','.join(matches)
    return result

## json 객체를 bed 형식으로
def jsontobed_byline(i):
    LocusId = i["LocusId"]
    LocusStructure = i["LocusStructure"]
    ReferenceRegion = i["ReferenceRegion"]
    # print(LocusId,LocusStructure,ReferenceRegion)
    # ReferenceRegion 여러개인 경우 list 처리가 됨
    # chrX:148500631-148500691"
    if type(ReferenceRegion) == list:
        #ReferenceRegion = '_'.join(ReferenceRegion)
        chrom, start = ReferenceRegion[0].split(':')
        start = start.split('-')[0]
        end = ReferenceRegion[-1].split('-')[-1]
    else:
        chrom, positions = ReferenceRegion.split(':')
        start, end = positions.split('-')
    # "LocusStructure": "(CAGG)*(CAGA)*(CA)*", to CAGG)n(CAGA)n(CA)n
    # "LocusStructure": "(GCC)*",
    STRUC = LocusStructure.replace("*","n").replace("+","n")
    MOTIFS= extract_motif(STRUC)
    bed_info = "ID=%s;MOTIFS=%s;STRUC=%s"%(LocusId,MOTIFS,STRUC)
    return "\t".join([chrom,start,end,bed_info])

def main():
    out = open(datain.replace(".json",".JSONtoBED.bed"),"w")
    for i in json_data:
        out.write(jsontobed_byline(i)+"\n")
    out.close()

main()


