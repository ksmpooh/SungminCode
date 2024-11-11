## STR annotation

library(tidyverse)
library(stringr)
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/db/annovar/")

anno_refseq <- read_table("output.vcf.avinput.variant_function",col_names = F) %>% select(1:4)
head(anno_refseq)
colnames(anno_refseq) <- c("type","gene","chrom","mid") 

anno_cyto <- read_table("output.vcf.avinput.hg38_cytoBand",col_names = F) %>% select(2:4) 

head(anno_cyto)

colnames(anno_cyto) <- c("cytoband","chrom","mid") 


anno_nested <- read_table("output.vcf.avinput.hg38_nestedRepeats",col_names = F) %>% select(2:4)
head(anno_nested)
colnames(anno_nested) <- c("nestedRepeat","chrom","mid") 
anno_nested %>% mutate(nestedRepeat = str_split_fixed(nestedRepeat,"=",2)[,2]) %>% 
  mutate(nestedRepeat = str_split_fixed(nestedRepeat,",",2)[,1]) %>% #count(nestedRepeat) %>% 
  mutate(nestedRepeat = ifelse(str_detect(nestedRepeat,"\\)n"),"simple",nestedRepeat)) %>% #count(nestedRepeat)
  mutate(nestedRepeat = ifelse(str_detect(nestedRepeat,"Alu"),"Alu",nestedRepeat)) %>% #count(nestedRepeat)
  mutate(nestedRepeat = ifelse(str_detect(nestedRepeat,"Arthur"),"Arthur",nestedRepeat)) %>% count(nestedRepeat)
mutate(nestedRepeat = ifelse(str_detect(nestedRepeat,"Arthur"),"Arthur",nestedRepeat)) %>% count(nestedRepeat)


ref <- read.table("output_bed_ref")
df_merge <- read.table("db.merge.all.bed")
df_merge %>% mutate(ID = str_replace_all(str_split_fixed(V4,";",3)[,1],"ID=","")) %>% select(-V4) -> df_merge
  

colnames(df_merge) <- c("chrom","start","end","STR_ID")
head(df_merge)
head(ref)

head()
head(anno_refseq)
head(anno_cyto)
ref %>% mutate(ID = paste0(V1,"-",V2)) -> ref
ref %>% select(4:5) %>% rename(STR_ID = V4) -> ref
head(ref)
ref %>% unique() %>% left_join(df_merge %>% unique()) -> ref
head(ref)

anno_refseq %>% mutate(ID = paste0(chrom,"-",mid)) -> anno_refseq
head(anno_refseq)

anno_cyto %>% mutate(ID = paste0(chrom,"-",mid)) -> anno_cyto
head(anno_cyto)
head(ref)

anno_refseq %>% select(type,gene,ID)%>% left_join(anno_cyto) 
anno_refseq %>% dim()

head(anno_refseq)
anno_refseq %>% unique() %>%  left_join(anno_cyto %>% unique()) %>% left_join(ref) -> anno
head(anno)
ref_cy <- read_table("hg38_cytoBand.txt",col_names = F)
head(ref_cy)
colnames(ref_cy) <- c("chrom","start","end","cytoband","cyto_type")
head(ref_cy)
ref_cy %>% filter(is.na(cyto_type))
table(ref_cy$cyto_type)
ref_cy %>% filter(cyto_type == 'acen') %>% mutate(arm = ifelse(str_detect(cytoband,"p"),"p","q")) %>% select(-cytoband,-cyto_type) -> ref_cy_centro
ref_cy %>% mutate(cytoband = paste0(str_replace_all(chrom,"chr",""),cytoband)) %>% #head()
  select(-start,-end) -> ref_cy
#  gpos 염색체 강도 -> heterochromatin
# gneg 유전자 밀도가 ???은 지역 ->, 전사가 잘 이러남 ,euchromatin
# gvar (Variable Region):
# acen (Centromeric Region):
# stalk : (acrocentric chromosomes): 주로 반복서열 , 13,14, 15 21 22 p-arm에 위치 satelite
head(ref_cy_centro)
head(ref_cy)
colnames(ref_cy_centro)<- c("chrom","acen_start","acen_end","arm")
ref_cy_centro %>% pivot_wider(names_from=arm,values_from = 2:3) %>% count(acen_start_p > acen_start_q)
ref_cy_centro %>% pivot_wider(names_from=arm,values_from = 2:3) %>% count(acen_end_p > acen_end_q)

ref_cy_centro %>% mutate(centro_pos_byarm = ifelse(arm == "p",acen_start,acen_end)) %>% select(-2,-3) -> ref_cy_centro
head(ref_cy_centro)
head(ref_cy)
head(anno)
anno %>% left_join(ref_cy) %>% head()
anno %>% left_join(ref_cy) %>% dim() # 321382     10

anno %>% left_join(ref_cy) %>% count(cyto_type)
anno %>% left_join(ref_cy) %>% mutate(arm = ifelse(str_detect(cytoband,"p"),"p","q")) %>% left_join(ref_cy_centro) %>% dim() #321382     10

anno %>% left_join(ref_cy)%>% mutate(arm = ifelse(str_detect(cytoband,"p"),"p","q")) %>% left_join(ref_cy_centro) -> anno1
anno1 %>% mutate(distance_bycen = ifelse(arm == "p",centro_pos_byarm - mid,mid-centro_pos_byarm)) %>% count(distance_bycen >= 0)
anno1 %>% count(type)


anno1 %>% #head()
  mutate(type = case_when(
    type == "upstream;downstream" ~ "upstream",
    type == "UTR5;UTR3" ~ "UTR5",
    type == "exonic;splicing" ~ "exonic",
    TRUE ~ type)) %>% #count(type)
  mutate(distance_bycen = ifelse(arm == "p",centro_pos_byarm - end,start-centro_pos_byarm)) %>%  #select(type,gene,ID,start,end)
  select(-ID,-arm,-centro_pos_byarm) %>% mutate(main_cyto_type = ifelse(str_detect(cyto_type,"gpos"),"gpos",cyto_type)) %>% #head()
  write.table("Raw.anno.merge.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")


anno1 %>% #head()
  mutate(type = case_when(
    type == "upstream;downstream" ~ "upstream",
    type == "UTR5;UTR3" ~ "UTR5",
    type == "exonic;splicing" ~ "exonic",
    TRUE ~ type)) %>% #count(type)
  mutate(distance_bycen = ifelse(arm == "p",centro_pos_byarm - end,start-centro_pos_byarm)) %>%  #select(type,gene,ID,start,end)
  select(-ID,-arm,-centro_pos_byarm) %>% mutate(main_cyto_type = ifelse(str_detect(cyto_type,"gpos"),"gpos",cyto_type)) %>%
  select(-gene,-chrom,-mid,-start,-end) %>%
  write.table("Raw.anno.merge.processing.onlyneed.txt",col.names = T,row.names = F,quote = F,sep = "\t")
