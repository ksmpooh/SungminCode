#ls *polished | xargs -I {} -P 80 bash -c "bash liftoff.after.sh {}"


input=$1
build=$2 #CHM13 or GRCh38

mkdir liftoff_processing
mkdir liftoff_processing/transcript/
mkdir liftoff_processing/gene/

awk '$3=="transcript" {
    attr = $NF  # 속성 필드는 보통 마지막 열에 있음
    source_gene = "NA"
    transcript_biotype = "NA"
    source_transcript = "NA"

    
    # 속성 필드를 ";"를 구분자로 배열에 저장
    n = split(attr, fields, ";")
    for (i = 1; i <= n; i++) {
        # 앞뒤 공백 제거 (선택사항)
        gsub(/^[ \t]+|[ \t]+$/, "", fields[i])
        if (fields[i] ~ /^gene_name=/) {
            split(fields[i], a, "=")
            gene_name = a[2]
        }
        else if (fields[i] ~ /^gene_biotype=/) {
            split(fields[i], a, "=")
            gene_biotype = a[2]
        }
        else if (fields[i] ~ /^source_transcript=/) {
            split(fields[i], a, "=")
            source_transcript = a[2]
        }
    }
    # $1은 보통 염색체 또는 시퀀스 ID
    print $1, gene_name, gene_biotype, source_transcript, source_transcript_name
}' $input > liftoff_processing/transcript/$input.afterliftoff_transcript





awk '$3=="gene" {
    attr = $NF  # 속성 필드는 보통 마지막 열에 있음
    gene_name = "NA"
    gene_biotype = "NA"
    extra_copy_number = "NA"
    coverage = "NA"
    
    # 속성 필드를 ";"를 구분자로 배열에 저장
    n = split(attr, fields, ";")
    for (i = 1; i <= n; i++) {
        # 앞뒤 공백 제거 (선택사항)
        gsub(/^[ \t]+|[ \t]+$/, "", fields[i])
        if (fields[i] ~ /^gene_name=/) {
            split(fields[i], a, "=")
            gene_name = a[2]
        }
        else if (fields[i] ~ /^gene_biotype=/) {
            split(fields[i], a, "=")
            gene_biotype = a[2]
        }
        else if (fields[i] ~ /^extra_copy_number=/) {
            split(fields[i], a, "=")
            extra_copy_number = a[2]
        }
        else if (fields[i] ~ /^coverage=/) {
            split(fields[i], a, "=")
            coverage = a[2]
        }
    }
    # $1은 보통 염색체 또는 시퀀스 ID
    print $1, gene_name, gene_biotype, extra_copy_number, coverage
}' $input > liftoff_processing/gene/$input.afterliftoff_gene







#1                   2       3       4       5       6   7   8
#KPPD001_h1tg000001l	Liftoff	gene	123176	141759	.	+	.	
#9
#ID=CHM13_G0003663;source_gene_common_name=FCGR2A;source_gene=ENSG00000143226.15;	
#10: gene_biotype
#protein_coding;gene_id=CHM13_G0003663;	
#11: gene_name
#FCGR2A;transcript_modes=transMap;Name=FCGR2A;source_transcript=N/A;alternative_source_transcripts=N/A;collapsed_gene_ids=N/A;collapsed_gene_names=N/A;paralogy=N/A;unfiltered_paralogy=N/A;alignment_id=N/A;frameshift=N/A;exon_anotation_support=N/A;intron_annotation_support=N/A;transcript_class=N/A;valid_start=N/A;valid_stop=N/A;proper_orf=N/A;extra_paralog=False;	
#12: coverage
#0.999;sequence_ID=0.998;valid_ORFs=4;	
#13: copy number
#0;copy_num_ID=CHM13_G0003663_0

########## coverage filter
## shell for geon_biotype
#mkdir gene_biotype

#gene : T2T -> GRCh38
#gene_name = "NA" ->  gene_name
#gene_biotype = "NA" -> gene_type
#extra_copy_number = "NA" -> extra_copy_number
#coverage = "NA" -> coverage


#T2T -> GRCh38
#source_gene -> Parent
#transcript_biotype = "NA" -> transcript_type
#source_transcript = "NA" -> transcript_id


    
#T2T
#source_transcript=ENST00000657469.1	source_transcript_name=LINC00708-203	source_gene=ENSG00000232170.6	
#GRCh38
#transcript_id=ENST00000778861.1; Parent=ENSG00000301427.1
#transcript_type

##T2T
### gene
#NA20129#2#JAHEPD010000001.1	Liftoff	gene	1428	10433	.	+	.	
#ID=CHM13_G0005600;source_gene_common_name=LINC00708;source_gene=ENSG00000232170.6;
#gene_biotype=lncRNA;gene_id=CHM13_G0005600;
#gene_name=LINC00708;transcript_modes=transMap;Name=LINC00708;
#source_transcript=N/A;alternative_source_transcripts=N/A;collapsed_gene_ids=N/A;collapsed_gene_names=N/A;paralogy=N/A;unfiltered_paralogy=N/A;alignment_id=N/A;frameshift=N/A;exon_anotation_support=N/A;intron_annotation_support=N/A;transcript_class=N/A;valid_start=N/A;valid_stop=N/A;proper_orf=N/A;extra_paralog=False;
#coverage=0.999;sequence_ID=0.998;
#extra_copy_number=0;copy_num_ID=CHM13_G0005600_0


##GRCh38
### gene
#HG01978#1#JAGYVS010000001.1	Liftoff	gene	1329	6545	.	-	.	
#ID=ENSG00000301427.1;
#gene_id=ENSG00000301427.1;
#gene_type=lncRNA;
#gene_name=ENSG00000301427;
#level=2;coverage=1.0;sequence_ID=0.995;extra_copy_number=0;copy_num_ID=ENSG00000301427.1_0

### transcript 
#HG01978#1#JAGYVS010000001.1	Liftoff	transcript	1329	6545	.	-	.	
#ID=ENST00000778861.1;Parent=ENSG00000301427.1;gene_id=ENSG00000301427.1;transcript_id=ENST00000778861.1;
#gene_type=lncRNA;gene_name=ENSG00000301427;
#transcript_type=lncRNA;transcript_name=ENST00000778861;
#level=2;tag=basic,Ensembl_canonical,TAGENE;
#extra_copy_number=0


#genome@genome205:/CDATA/pangenome/HPRC/annotation/gencode_CHM13$ awk '$3=="gene"{print $9}' NA20129.2.gencode_CHM13.liftoff.gff3_polished | sed 's/;/\t/g' |head
#ID=CHM13_G0005600	 1
#source_gene_common_name=LINC00708	2
#source_gene=ENSG00000232170.6	3
#gene_biotype=lncRNA	4
#gene_id=CHM13_G0005600	5
#gene_name=LINC00708	6
#transcript_modes=transMap 7
#Name=LINC00708 8
#source_transcript=N/A 9
#alternative_source_transcripts=N/A	collapsed_gene_ids=N/A	10
#collapsed_gene_names=N/A	11
#paralogy=N/A 12
#unfiltered_paralogy=N/A 13
#alignment_id=N/A 14
#frameshift=N/A 15
#exon_anotation_support=N/A
#intron_annotation_support=N/A	transcript_class=N/A	valid_start=N/A	valid_stop=N/A	proper_orf=N/A	extra_paralog=False	coverage=0.999	sequence_ID=0.998	extra_copy_number=0	copy_num_ID=CHM13_G0005600_0

#genome@genome205:/CDATA/pangenome/HPRC/annotation/gencode_CHM13$ awk '$3=="transcript"{print $9}' NA20129.2.gencode_CHM13.liftoff.gff3_polished | sed 's/;/\t/g' | head -n 1
#ID=CHM13_T0020871	source_transcript=ENST00000657469.1	source_transcript_name=LINC00708-203	source_gene=ENSG00000232170.6	
#transcript_modes=transMap	gene_biotype=lncRNA	transcript_biotype=lncRNA	alignment_id=ENST00000657469.1-0	frameshift=nan	
#exon_annotation_support=1,1,1,1,1	intron_annotation_support=1,1,1,1	transcript_class=ortholog	valid_start=True	valid_stop=True	adj_start=nan	adj_stop=nan	
#proper_orf=True	level=2	hgnc_id=HGNC:44694	tag=basic,TAGENE	
#havana_gene=OTTHUMG00000017645.2	havana_transcript=OTTHUMT00000516776.1	paralogy=nan	unfiltered_paralogy=nan	source_gene_common_name=LINC00708	transcript_id=CHM13_T0020871	
#gene_id=CHM13_G0005600	Parent=CHM13_G0005600	
#transcript_name=LINC00708-203	Name=LINC00708	gene_name=LINC00708	
#alternative_source_transcripts=N/A	collapsed_gene_ids=N/A	collapsed_gene_names=N/A	extra_paralog=False	extra_copy_number=0