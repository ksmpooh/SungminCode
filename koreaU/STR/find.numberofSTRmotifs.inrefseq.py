'''
DAB1
NOTCH2NLA
NOTCH2NLC
NUTM2B-AS1
FRA10AC1
C11ORF80
ATXN8
ZIC2
XYLT1
TNRC6A
BEAN1
COMP
STARD7
AFF3
HOXD13
FOXL2
YEATS2
RAPGEF2
MARCHF6
RUNX2
HOXA13
ZNF713
LRP12
SAMD12
ARX
SOX3
TMEM185A
'''
import os,glob

fileIn = "pathogenic_repeats.hg38.bed"
df = open(fileIn,"r")
df = [s.strip() for s in df]

target_ID = ["DAB1","NOTCH2NLA","NOTCH2NLC","NUTM2B-AS1","FRA10AC1","C11ORF80","ATXN8","ZIC2","XYLT1","TNRC6A","BEAN1","COMP","STARD7","AFF3","HOXD13","FOXL2","YEATS2","RAPGEF2","MARCHF6","RUNX2","HOXA13","ZNF713","LRP12","SAMD12","ARX","SOX3","TMEM185A"]

for i in df:
    tmp = i.split("\t")
    chrom = tmp[0]
    start = tmp[1]
    end = tmp[2]
    id = tmp[3].split(";")[0].split('=')[1]
    if id not in target_ID:
        continue
    print("%s:%s-%s\t%s"%(chrom,start,end,id))



for i in df:
    tmp = i.split("\t")
    chrom = tmp[0]
    start = tmp[1]
    end = tmp[2]
    id = tmp[3].split(";")[0].split('=')[1]
    motifs = tmp[3].split(";")[1]
    motif = motifs.split('=')[1]
    if id not in target_ID:
        continue
    if "," in motifs:
        print("%s:%s-%s\t%s\t%s"%(chrom,start,end,id,motif))



#/BDATA/smkim/HLAseq/REF/human_GRCh38_no_alt_analysis_set
    
#awk '{print "samtools faidx ../human_GRCh38_no_alt_analysis_set.basic.fasta "$1" > refseq/"$2"_"$3".txt"}' onlytrgt.pathoID.txt
#RFC1 G/A -> R
'''
    {
        "LocusId": "ATXN8OS",
        "LocusStructure": "(CTA)*(CTG)*",
        "ReferenceRegion": [
            "chr13:70139353-70139383",
            "chr13:70139383-70139428"
            chr13:70139353-70139428
>chr13:70139353-70139428
A CTACTACTACTACTACTACTACTACTACTA CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG

'''


'''
########## good case (예시, trgt에도 있고, eh에도 있고)
trgt bed: chr13	70139353	70139428	ID=ATXN8;MOTIFS=CTA,CTG;STRUC=(CTA)n(CTG)n
eh json {
        "LocusId": "ATXN8OS",
        "LocusStructure": "(CTA)*(CTG)*",
        "ReferenceRegion": [
            "chr13:70139353-70139383",
            "chr13:70139383-70139428"}

                
#ref region 
>chr13:70139353-70139428
A CTACTACTACTACTACTACTACTACTACTA CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
1 30                             45
70139353 70139353+30=70139383    70139383+45 = 70139428

########### bed case (예시, trgt에 O, eh에 X)
trgt bed: chr16	66490398	66490467	ID=BEAN1;MOTIFS=TGGAA,TAAAA;STRUC=(TGGAA)n(TAAAA)n

#ref region 
>chr16:66490398-66490467
A TAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA TAA TAAAATAAAAA
1 55                                                      3   11
????????


'''
==> ATXN8_CTA,CTG.txt <==
>chr13:70139353-70139428
A CTACTACTACTACTACTACTACTACTACTA CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG

==> BEAN1_TGGAA,TAAAA.txt <== # 추가 확인 필요 TRRAA -> (TGGAA)*(TAAAA)*
>chr16:66490398-66490467
ATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAATAAAATAAAAA

==> DAB1_AAAAT,GAAAT.txt <==
>chr1:57367043-57367119
C AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA

==> MARCHF6_TTTTA,TTTCA.txt <==
>chr5:10356346-10356412
TATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTC

==> RAPGEF2_TTTTA,TTTCA.txt <==
>chr4:159342526-159342617
CTTTTATTTTATTTTATTTTATTTTATATTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTAT

==> SAMD12_TGAAA,TAAAA.txt <== 
>chr8:118366815-118366919
ATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAATAAAATAAAATAAAATAAAATAAAATAAAATAAAAA

==> STARD7_TGAAA,TAAAA.txt <==
>chr2:96197066-96197122
CAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATA

==> TNRC6A_TTTTA,TTTCA.txt <==
>chr16:24613439-24613530
ATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTATTTTAT

==> YEATS2_TTTTA,TTTCA.txt <==
>chr3:183712187-183712223
CTTTTATTTTATTTTATTTTATTTTATTTTATTTTAT