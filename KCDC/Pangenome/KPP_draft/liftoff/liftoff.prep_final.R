#### lifteroff finael
#cd /ADATA/pangenome/liftoff/with_flagger
library(tidyverse)

# 106

multi_exon_gene_list <- read_table("/CDATA/pangenome/liftoff/db/gencode/pc_multi_exon_gene_names.txt",col_names = F)
df <- read_table("/ADATA/pangenome/liftoff/with_flagger/merge/merge.with_flagger.gene_protein_coding.Hap.selected.txt") %>%
  mutate(ID = str_split_fixed(contig,"tg",2)[,1])
df %>% filter(Flagger == "Hap") %>% count(tag)

df %>% filter(str_detect(Flagger,"Hap")) %>% count(tag)


df_all <- read_table("/ADATA/pangenome/liftoff/with_flagger/all_merge/merge.with_flagger.gene_protein_coding.selected.txt") %>%
  mutate(ID = str_split_fixed(contig,"tg",2)[,1])

df_all_withhgnc <- read_table("/ADATA/pangenome/liftoff/with_flagger/all_merge/merge.with_flagger.gene_protein_coding.selected_withhgnc.txt") %>%
  mutate(ID = str_split_fixed(contig,"tg",2)[,1])

df_all_withhgnc %>% filter(!is.na(hgnc_id)) %>% 
  #filter(Flagger == "Hap") %>% 
  count(ID,hgnc_id) %>% filter(n != 1) %>%  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

df_all_withhgnc %>% filter(!is.na(hgnc_id)) %>% filter(Flagger == "Hap") %>% count(ID,hgnc_id) %>% filter(n != 1) %>%  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
  #filter(is.na(tag)) %>% 
  filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  #filter(valid_ORFs > 0) %>% 
  #filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
  


df_all %>% select(ID,gene_name,extra_copy_number) %>%
  group_by(ID,gene_name) %>%
  summarise(
    extra_copy_number = max(extra_copy_number, na.rm = TRUE),
    .groups = "drop") -> df_all_drop
  
df_all_drop %>% filter(extra_copy_number != 0) %>% count(ID)
df_all_drop %>% filter(extra_copy_number != 0) %>% count(gene_name) %>% dim
df_all_drop %>%  filter(extra_copy_number != 0) %>% filter(gene_name == "RFPL4AL1")
df_all_drop %>% filter(extra_copy_number == 3) 
df_all %>% select(ID,gene_name,extra_copy_number) %>% filter(ID == "KPPD001_h1", gene_name == "RFPL4AL1")
head(df) 
df %>% filter(str_detect(Flagger,"Hap"),
              is.na(tag) | tag == "confirm_experimentally") %>%  count(Flagger,tag)


df %>% filter(str_detect(Flagger,"Hap"),
              is.na(tag) | tag == "confirm_experimentally") %>%  count(coverage > 0.9,sequence_ID > 0.9)


df %>% filter(str_detect(Flagger,"Hap"),
              is.na(tag) | tag == "confirm_experimentally") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% count(level)


df %>% filter(str_detect(Flagger,"Hap"),
              is.na(tag) | tag == "confirm_experimentally") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% count(Flagger)


df %>% filter(Flagger == "Hap" | str_detect(Flagger,"Hap,"),
              is.na(tag) | tag == "confirm_experimentally") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% filter(str_detect(Flagger,",")


df %>% filter(str_detect(Flagger,"Hap"),
              is.na(tag) | tag == "confirm_experimentally") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  filter(extra_copy_number != 0) %>%
  mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% count(ID,gene_name) %>% count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
`mean(n)` `median(n)` `min(n)` `max(n)`
<dbl>       <dbl>    <int>    <int>
  1      56.8        56.5       28       84
'''

df %>%  filter(extra_copy_number != 0) %>% 
  filter(Flagger == "Hap",is.na(tag) | tag == "confirm_experimentally") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% count(ID,gene_name) %>% count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
'''
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      54.2          53       28       80
'''

df %>%  filter(extra_copy_number != 0) %>% 
  filter(Flagger == "Hap",is.na(tag) | tag == "confirm_experimentally") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% count(ID,gene_name) %>% count(ID) %>% summarise(mean(n),median(n),min(n),max(n))


df %>% filter(Flagger == "Hap",is.na(tag) | tag == "confirm_experimentally") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% 
  filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% 
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
`mean(n)` `median(n)` `min(n)` `max(n)`
<dbl>       <dbl>    <int>    <int>
  1      33.0          32       17       53
'''
#0	1	1	1	1	57.750000	56.000000	41	87	132

df %>%  
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% 
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))


df %>%  
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  mutate(ID = str_split_fixed(contig,"tg",2)[,1]) -> df_QC

df %>%  
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      71.8          72       46      116
'''

df %>%  filter(gene_name  %in% multi_exon_gene_list$X1) %>%
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      49.9          49       28       80
'''
  
df %>% filter(gene_name  %in% multi_exon_gene_list$X1) %>%
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
# A tibble: 1 × 4
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      44.0          43       23       71
'''
  
df %>% filter(gene_name  %in% multi_exon_gene_list$X1) %>%
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  count(ID,gene_name) %>% filter(n != 1) %>% count(gene_name)
# 947

df %>% filter(gene_name  %in% multi_exon_gene_list$X1) %>%
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
#  filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  count(ID,gene_name) %>% filter(n != 1) %>% count(gene_name)
#989
  

df %>% filter(gene_name  %in% multi_exon_gene_list$X1) %>%
  filter(Flagger == "Hap") %>%  filter(coverage >= 0.9,sequence_ID >= 0.9) %>% 
  filter(valid_ORFs > 0) %>% filter(level < 3) -> df_QC


df_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% write.table("/ADATA/pangenome/liftoff/with_flagger/summary/liftoff.CNV.bysample.txt",col.names = T,row.names = F,quote = F,sep = "\t")

df_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(gene_name) %>% write.table("/ADATA/pangenome/liftoff/with_flagger/summary/liftoff.CNV.byGene.txt",col.names = T,row.names = F,quote = F,sep = "\t")

gene_fam <- read.table("/CDATA/pangenome/liftoff/db/gencode/gencode_pc_gene_to_hgnc_family.primary.tsv",sep="\t",header=T)


gene_id_mat <- df_QC %>%
  count(ID, gene_name) %>%
  pivot_wider(
    names_from  = ID,
    values_from = n,
    values_fill = 0
  )

df_QC %>%
  count(ID, gene_name) %>% filter(n != 1) %>% 
  mutate(n = n-1) %>%
  pivot_wider(names_from  = ID,values_from = n,values_fill = 0) %>%
  arrange(gene_name) %>%   # ← 여기
  write.table("/ADATA/pangenome/liftoff/with_flagger/summary/liftoff.CNV.raw.txt",col.names = T,row.names = F,quote = F,sep = "\t")



df_QC_fam <- df_QC %>%
  left_join(gene_fam %>% 
              select(gene_name,primary_family_id, primary_family_abbreviation, primary_family_name),by = "gene_name"
  )

fam_summary <- df_QC_fam %>%
  filter(!is.na(primary_family_id)) %>%
  group_by(primary_family_id, primary_family_abbreviation, primary_family_name) %>%
  summarise(
    n_records = n(),                         # 전체 행 수 (샘플/contig 단위 포함)
    n_genes   = n_distinct(gene_name),       # family 내 고유 gene 수
    n_contigs = n_distinct(contig),          # 관측된 contig 수
    n_flagger_classes = n_distinct(Flagger), # Hap/Dup/Col/Err 종류 수
    mean_cov = mean(coverage, na.rm = TRUE),
    median_cov = median(coverage, na.rm = TRUE),
    mean_seqID = mean(sequence_ID, na.rm = TRUE),
    mean_validORF = mean(valid_ORFs, na.rm = TRUE),
    mean_extraCN  = mean(extra_copy_number, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genes), primary_family_id)
  #"/ADATA/pangenome/liftoff/with_flagger/summary"
  
fam_summary
write.table(fam_summary,"/ADATA/pangenome/liftoff/with_flagger/summary/liftoff.CNV.GeneFam.sum.txt",col.names = T,row.names = F,quote = F,sep = "\t")

# 53개(범위 27–100) gene = 1,367

cpc <- read_table("/CDATA/pangenome/liftoff/merge/CPC/CPC.all.merge.with_flagger.gene_protein_coding.select.txt") %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])
cpc_genefam <- read_table("/CDATA/pangenome/liftoff/merge/CPC/CPC.all.merge.with_flagger.gene_protein_coding.selected_withhgnc.txt") %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])
cpc %>% filter(is.na(tag)) %>% count(coverage >= 0.9,sequence_ID>=0.9)

cpc_genefam %>% count(ID,hgnc_id) %>% filter(n != 1) %>% count(ID)  %>% 
  summarise(mean(n),median(n),min(n),max(n))



cpc_genefam %>%
  #filter(is.na(tag)) %>% 
  #filter(coverage >= 0.9) %>%
  #filter(sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

cpc %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% 
  filter(gene_name %in% multi_exon_gene_list$X1) %>%
  #filter(is.na(tag)) %>% 
  filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  #filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      42.4          42       20       89
'''

cpc %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% 
  filter(is.na(tag)) %>% 
  filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% #filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      43.2          43       21       90
'''
cpc %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% 
  filter(is.na(tag)) %>% 
  filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  #filter(valid_ORFs > 0) %>% #filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
'''

# A tibble: 1 × 4
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      58.1        59.5       31      113
'''
cpc_gene_paper <- read_table("/CDATA/pangenome/liftoff/db/cpc.gene.list.txt")
cpc %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% 
  filter(gene_name %in% cpc_gene_paper$gene) %>% #count(gene_name) %>% arrange(-n)
  #filter(is.na(tag)) %>% 
  #filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  #filter(valid_ORFs > 0) %>% 
  #filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))


'''
# A tibble: 1 × 4
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      66.5          67       37      118
'''



cpc %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% 
  filter(gene_name  %in% multi_exon$gene_name) %>%
  #filter(is.na(tag)) %>% 
  filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
'''
'''

cpc %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1]) %>% #head
  filter(extra_copy_number == 1) %>%
  #filter(is.na(tag)) %>% 
  #filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  #filter(valid_ORFs > 0) %>% filter(level < 3) %>%
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
'''
# A tibble: 1 × 4
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      81.3        81.5       48      143
'''

cpc %>% 
  filter(is.na(tag)) %>% 
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  #filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))


cpc %>% filter(extra_copy_number == 1) %>% count(ID,gene_name) %>% count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

cpc %>% filter(extra_copy_number > 0) %>% 
  count(ID,gene_name) %>% count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

##### new DB 20260112
library(tidyverse)
df_new <- read_table("/ADATA/pangenome/liftoff/01.liftoff/gencode_GRCh38.v38_partial/gene/with_flagger/merge.liftoff.gene.withflagger.txt") 
df_new <- df_new %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])


df_new %>% count(coverage >= 0.9)
df_new %>% count(valid_ORFs > 0)

df_gene <- read_table("/ADATA/pangenome/liftoff/02.annotation/merge.genecode.v38.prep.txt_forR") %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])
df_trans <- read_table("/ADATA/pangenome/liftoff/02.annotation/merge.genecode_transcript.v38.prep.txt_forR") %>% mutate(ID = str_split_fixed(contig,"tg",2)[,1])

grch_gene <- read_table("/BDATA/smkim/pangenome/liftoff/db/gencode.v38.primary_assembly.annotation.gff3_gene_prep")
grch_trans <- read_table("/BDATA/smkim/pangenome/liftoff/db/gencode.v38.primary_assembly.annotation.gff3_trans")


grch_gene %>% count(gene_type == "protein_coding")
'''
`gene_type == "protein_coding"`     n
<lgl>                           <int>
  1 FALSE                           40724 -> 39482
2 TRUE                            19984 -> 39483
'''
grch_trans %>% count(transcript_type == "protein_coding")


'''
  `transcript_type == "protein_coding"`      n
  <lgl>                                  <int>
1 FALSE                                 150285
2 TRUE                                   86794
'''
df_trans %>% mutate(coding_type = ifelse(transcript_type=="protein_coding","Protein_coding_trans","Noncoding_trans")) %>%
  select(ID,transcript_name,coding_type) %>% unique() %>%
  count(ID,coding_type) %>% write.table("/ADATA/pangenome/liftoff/99.summary/transcript.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 


df_gene %>% mutate(coding_type = ifelse(gene_type=="protein_coding","Protein_coding_genes","Noncoding_genes")) %>% 
  select(ID,gene_name,coding_type) %>% unique() %>%
  count(ID,coding_type) %>% write.table("/ADATA/pangenome/liftoff/99.summary/gene.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 



df_new %>% filter(extra_copy_number != 0) %>% count(gene_type)

df_new %>% filter(gene_type == "protein_coding") %>% count(valid_ORFs)
df_new %>% filter(gene_type == "protein_coding") %>% count(valid_ORFs != 0)
df_new %>% filter(gene_type == "protein_coding") %>% filter(is.na(valid_ORFs))

df_new %>% filter(gene_type == "protein_coding") -> df_new_proteincoding

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
\
'''
# A tibble: 1 × 4
`mean(n)` `median(n)` `min(n)` `max(n)`
<dbl>       <dbl>    <int>    <int>
  1      67.8          67       43      111
'''



df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  #filter(valid_ORFs > 0) %>% 
  #filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))



df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>% count(ID,gene_name) %>%  filter(n != 1) %>% 
  count(gene_name) %>% dim()
#[1] 1080    2

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

'''
`mean(n)` `median(n)` `min(n)` `max(n)`
<dbl>       <dbl>    <int>    <int>
  1      49.8          49       30       83
'''

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(gene_name  %in% multi_exon$gene_name) %>%
  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>%filter(valid_ORFs > 0) %>%  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
'''
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      40.9          40       22       69
'''
df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>% #filter(valid_ORFs > 0) %>%  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
'''
# A tibble: 1 × 4
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      67.1          67       42      110
'''

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(gene_name) %>% dim()
# 1154    2

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>% #filter(valid_ORFs > 0) %>%  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(gene_name) %>% dim()
#[1] 1152    2


df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>% filter(valid_ORFs > 0) %>%  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(gene_name) %>% dim()
# 1080    2

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>% filter(valid_ORFs > 0) %>%  filter(!is.na(valid_ORFs)) %>%
  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(gene_name) %>% dim()

#[1] 1094    2

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>% 
  filter(valid_ORFs > 0) %>%  #filter(!is.na(valid_ORFs)) %>%
  #filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(gene_name) %>% dim()


df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(gene_name  %in% multi_exon$gene_name) %>%
  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>% filter(valid_ORFs > 0) %>%  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(gene_name) %>% dim()
#940

df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(gene_name  %in% multi_exon$gene_name) %>%
#  filter(coverage >= 0.9) %>% filter(sequence_ID>=0.9) %>% filter(valid_ORFs > 0) %>%  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(gene_name) %>% dim()
#980   2
  
df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>%
  count(ID,gene_name) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))
'''
  `mean(n)` `median(n)` `min(n)` `max(n)`
      <dbl>       <dbl>    <int>    <int>
1      49.8          49       30       83
'''
df_new_proteincoding %>% filter(Flagger == "Hap") %>%
  filter(coverage >= 0.9) %>%
  filter(sequence_ID>=0.9) %>%
  filter(level < 3) %>%
  filter(valid_ORFs > 0) -> df_new_proteincoding_QC
'''
`mean(n)` `median(n)` `min(n)` `max(n)`
<dbl>       <dbl>    <int>    <int>
  1      6.27           6        4       12
  '''

df_new_proteincoding %>% 
  count(ID,hgnc_id) %>%  filter(n != 1) %>% #head
  count(ID) %>% summarise(mean(n),median(n),min(n),max(n))

df_new_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% write.table("/ADATA/pangenome/liftoff/99.summary/liftoff.CNV.bysample.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table("/ADATA/pangenome/liftoff/99.summary/gene.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 
df_new_proteincoding_QC %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(gene_name) %>% write.table("/ADATA/pangenome/liftoff/99.summary/liftoff.CNV.byGene.txt",col.names = T,row.names = F,quote = F,sep = "\t")

df_new_proteincoding_QC %>%
  count(ID, gene_name) %>% filter(n != 1) %>% 
  mutate(n = n-1) %>%
  pivot_wider(names_from  = ID,values_from = n,values_fill = 0) %>%
  arrange(gene_name) %>%   # ← 여기
  write.table("/ADATA/pangenome/liftoff/99.summary/liftoff.CNV.raw.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df_new_proteincoding_QC %>% filter(gene_name  %in% multi_exon$gene_name) %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(ID) %>% write.table("/ADATA/pangenome/liftoff/99.summary/liftoff.CNV_multiexon.bysample.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table("/ADATA/pangenome/liftoff/99.summary/gene.annotation.byhap.txt",col.names = T,row.names = F,quote = F,sep = "\t") 
df_new_proteincoding_QC %>% filter(gene_name  %in% multi_exon$gene_name) %>%
  count(ID,gene_name) %>% filter(n != 1) %>% 
  count(gene_name) %>% write.table("/ADATA/pangenome/liftoff/99.summary/liftoff.CNV_multiexon.byGene.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df_new_proteincoding_QC %>% filter(gene_name  %in% multi_exon$gene_name) %>%
  count(ID, gene_name) %>% filter(n != 1) %>% 
  mutate(n = n-1) %>%
  pivot_wider(names_from  = ID,values_from = n,values_fill = 0) %>%
  arrange(gene_name) %>%   # ← 여기
  write.table("/ADATA/pangenome/liftoff/99.summary/liftoff.CNV_multiexon.raw.txt",col.names = T,row.names = F,quote = F,sep = "\t")



multi_exon <- read_table("/BDATA/smkim/pangenome/liftoff/db/genocode.v38.multiexon_genes.tsv") %>% filter()
hprc_gene_list <- read_table("/CDATA/pangenome/liftoff/db/hprc.inDB.info.txt")


grch <- read_table("/CDATA/pangenome/liftoff/db/gencode/ref_liftoff/gencode_GRCh38.v38_partial/GRCh38.primary_assembly.genome.gencode_GRCh38.v38_partial.liftoff.gff3_polished_gene.prep")

grch %>% filter(gene_type == "protein_coding") %>%
  filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>% filter(contig == "KI270731.1")

grch %>% filter(gene_type == "protein_coding") %>%
  filter(coverage >= 0.9,sequence_ID>=0.9) %>%
  filter(valid_ORFs > 0) %>% 
  filter(level < 3) %>% #filter(gene_name == "GGTLC3")count(contig) -> a
  count(gene_name) -> grch_cnv

hprc_gene_list %>% select(gene,GRCh38) -> hprc_grch_38
grch_cnv %>% merge(hprc_grch_38,by.x = "gene_name",by.y = "gene")
grch_cnv %>% merge(hprc_grch_38,by.x = "gene_name",by.y = "gene") %>% filter(n-1 != GRCh38) %>% head
