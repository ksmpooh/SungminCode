## pathogenic 20240826
library(tidyverse)
library(data.table)
library(ggbreak)
library(ggpubr)
library(ggplot2)
library(cowplot)

## server에서 작업 후

load(file="~/Desktop/KU/@research/STR/02.compare/STR_prep/eh.ganstr.v13.raw.merge.RDATA")
load(file="~/Desktop/KU/@research/STR/02.compare/STR_prep/trgt.ganstr.v13.raw.merge.RData")

###ref
dup_str <- read.table("~/Desktop/KU/@research/STR/figure/rm_dup_str_list.txt",header = T)
head(dup_str)
ru <- read_table("~/Desktop/KU/@research/STR/db/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed",col_names = F)
head(ru)
ru %>% mutate(MOTIFS = str_split_fixed(X4,";",3)[,2]) %>% 
  mutate(MOTIFS = str_split_fixed(MOTIFS,"=",2)[,2]) %>% #head()
  mutate(ID = str_split_fixed(X4,";",3)[,1]) %>% 
  mutate(ID = str_split_fixed(ID,"=",2)[,2]) %>% #head()
  select(X1,X2,X3,MOTIFS,ID) -> ru
head(ru)

ref <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/db/str-analysis/str_analysis/variant_catalogs/catalog.GRCh38.with_adjacent_repeats.TRGT.bed",col_names = F)
head(ref)
ref %>%  mutate(MOTIFS = str_split_fixed(X4,";",3)[,2]) %>% 
  mutate(MOTIFS = str_split_fixed(MOTIFS,"=",2)[,2]) %>% #head()
  mutate(ID = str_split_fixed(X4,";",3)[,1]) %>% 
  mutate(ID = str_split_fixed(ID,"=",2)[,2]) %>% #head()
  select(X1,X2,X3,MOTIFS,ID) -> ref
head(ref)
ref$ID2 <- paste0(ref$X1,"_",ref$X2,"_",ref$X3)
ru$ID2 <- paste0(ru$X1,"_",ru$X2,"_",ru$X3)

head(ru)
ru %>% filter(str_detect(ID,"chr"))%>% 
  filter(ID2 %in% ref$ID2) -> rmlist2

table(ref$ID2 %in% ru$ID2)
head(ru)
head(ref)
dup_str %>% select(ID) %>% rbind(rmlist2  %>% select(ID)) -> rmlist_final
head(rmlist_final)
head(dup_str)
head(ru);dim(ru)
head(ref);dim(ref)
ru %>% filter(str_detect(ID,"chr")) %>% filter(!(ID %in%rmlist_final$ID)) %>% dim() #321253
321253 + 73 #321326
ru %>% filter(str_detect(ID,"chr")) %>% filter(!(ID %in%rmlist_final$ID)) %>% rbind(ref) %>% dim() #321326 #73
ru %>% filter(str_detect(ID,"chr")) %>% filter(!(ID %in%rmlist_final$ID)) %>% rbind(ref) %>% select(ID2) %>% unique() %>% dim #321326
ru %>% filter(str_detect(ID,"chr")) %>% filter(!(ID %in%rmlist_final$ID)) %>% rbind(ref) -> final_ref
colnames(final_ref) <- c("chrom","start","end","MOTIFS","ID","STR_region")
write.table(final_ref,"~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",col.names = T,row.names = F,quote = F,sep = "\t")
###
head(eh)
head(trgt)
head(eh_patho)
head(trgt_patho)

eh_patho <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_prep/eh.pathogenic_analysisDB.raw.merge.txt",header = T)
trgt_patho <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_prep/trgt.pathogenic_analysisDB.raw.merge.txt",header = T)
#NIH23J3904558	NA	NIH20N2594890

eh %>% select(ID,STR_ID,FILTER,type) %>% filter(ID != "NIH20N2594890") %>% filter(str_detect(STR_ID,"chr"))-> eh_qc
trgt %>% select(ID,TRID,Allele,MC,AP) %>% filter(ID != "NIH23J3904558") %>% filter(str_detect(TRID,"chr"))-> trgt_qc
head(trgt_qc)
eh_patho %>% select(ID,STR_ID,FILTER,type) %>% filter(ID != "NIH20N2594890") %>% filter(str_detect(STR_ID,"chr"))
head(eh_patho)
trgt_qc %>% count(MC == "0")
trgt_qc %>% count(AP == ".")
trgt %>% filter(is.na(MC))
'''
   MC == "0"        n
      <lgcl>    <int>
1:     FALSE 41576918
2:      TRUE    24452
3:        NA   164770
> trgt_qc %>% count(AP == ".")
   AP == "."        n
      <lgcl>    <int>
1:     FALSE 41586447
2:      TRUE    14923
3:        NA   164770
'''
trgt_qc %>% filter(AP != '.') %>% count(TRID) -> trgt_qc_count
eh_qc %>% filter(FILTER == "PASS") %>% count(STR_ID) -> eh_qc_count

head(eh_qc)
eh_qc %>% count(FILTER)

VennDiagram::venn.diagram(x = list(long = trgt_qc_count %>% filter(n == 130) %>% rownames(),
                                   short = eh_qc_count %>% filter(n == 65) %>% rownames()), filename = NULL,
             alpha = 0.3, # 원의 반투명도를 설정
        고리 레이블 크기
             cat.cex = 1.9, # 집합 이름의 글자 크기
             cat.fontface = "bold", # 글자 폰트 스타일 설정
             cat.fontfamily = "Arial") -> p
cowplot::plot_grid(p)
head(trgt_qc)
trgt_qc %>% filter(AP != '.') %>% count(TRID) -> trgt_qc_count
eh_qc %>% filter(FILTER == "PASS") %>% filter(type == "SPANNING/SPANNING")%>% count(STR_ID) -> eh_qc_count

eh_qc %>% filter(FILTER == "PASS") %>% count(type)
head(eh_qc_count)
head(eh_qc)
eh_qc %>% count(FILTER)

VennDiagram::venn.diagram(x = list(long = trgt_qc_count %>% filter(n == 130) %>% rownames(),
                                   short = eh_qc_count %>% filter(n == 65) %>% rownames()), filename = NULL,
                          #category.names = c("STR DB", "Long-read", "Short-read"),
                          #col = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원의 색상을 반투명하게 설정
                          #fill = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원 내부 색상도 반투명하게 설정
                          alpha = 0.3, # 원의 반투명도를 설정
                          cex = 1.9, # 카테고리 레이블 크기
                          cat.cex = 1.9, # 집합 이름의 글자 크기
                          cat.fontface = "bold", # 글자 폰트 스타일 설정
                          cat.fontfamily = "Arial") -> p
cowplot::plot_grid(p)

head(trgt_patho)
head(eh_patho)
eh_patho %>% select(ID,STR_ID,FILTER,TYPE) %>% filter(ID != "NIH20N2594890") -> eh_patho_qc
trgt_patho %>% select(ID,TRID,Allele,MC,AP) %>% filter(ID != "NIH23J3904558") -> trgt_patho_qc


#eh_patho %>% select(ID,STR_ID,FILTER) -> eh_patho_qc
#trgt_patho %>% select(ID,TRID,Allele,MC,AP) -> trgt_patho_qc
head(trgt_qc)
trgt_patho_qc %>% count(MC == "0")
trgt_patho_qc %>% count(AP == ".")
trgt_patho_qc %>% filter(is.na(MC))
trgt_patho %>% filter(is.na(MC))
'''
  MC == "0"    n
1     FALSE 9500
2      TRUE  114
3        NA   22
> trgt_patho_qc %>% count(AP == ".")
  AP == "."    n
1     FALSE 9614
2        NA   22
'''
eh_patho_qc %>% count(FILTER)
trgt_patho_qc %>% filter(AP != '.') %>% count(TRID) -> trgt_patho_qc_count
eh_patho_qc %>% filter(FILTER == "PASS") %>% count(STR_ID) -> eh_patho_qc_count

VennDiagram::venn.diagram(x = list(long = trgt_patho_qc_count %>% filter(n == 132) %>% rownames(),
                                   short = eh_patho_qc_count %>% filter(n == 66) %>% rownames()), filename = NULL,
                          #category.names = c("STR DB", "Long-read", "Short-read"),
                          #col = c("#440154AA", "#21908DAA", "#F하게 설정
                          #fill = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원 내부 색상도 반투명하게 설정
                          alpha = 0.3, # 원의 반투명도를 설정
                          cex = 1.9, # 카테고리 레이블 크기
                          cat.cex = 1.9, # 집합 이름의 글자 크기
                          cat.fontface = "bold", # 글자 폰트 스타일 설정
                          cat.fontfamily = "Arial") -> p
cowplot::plot_grid(p)
head(short)


#####

head(eh_qc)
head(trgt_qc)
head(trgt_patho_qc)
head(eh_patho_qc)

eh_qc$type <- "simple"
trgt_qc$type <- "simple"

eh_qc %>% rbind(eh_patho_qc) -> eh_merge_qc
trgt_qc %>% rbind(trgt_patho_qc) -> trgt_merge_qc

head(eh_qc)
eh_qc %>% filter(STR_ID %in% rmlist_final$ID)
trgt_merge_qc
eh_merge_qc
#dim(table(eh_merge_qc$ID))
eh_merge_qc %>% group_by(type,ID) %>% summarise(n=prop.table(FILTER))
#trgt_patho_qc %>% 

  
head(eh_qc)
head(trgt_qc)
eh_qc %>% filter(!(STR_ID %in% rmlist_final$ID)) -> eh_qc_rmlist
trgt_qc %>% filter(!(TRID %in% rmlist_final$ID)) -> trgt_qc_rmlist

trgt_patho_qc
eh_qc_rmlist
head(ref)
trgt_qc_rmlist
trgt_qc_rmlist %>% filter(is.na(MC))
trgt_qc_rmlist %>% count(TRID %in% ref$ID2)
eh_qc_rmlist %>% count(STR_ID %in% ref$ID2) 

eh_qc_rmlist %>% filter(FILTER == "PASS") %>% count(STR_ID) %>% filter(n==65) %>% dim #314171
eh_qc_rmlist %>% filter(FILTER == "PASS") %>% filter(type == "SPANNING/SPANNING") %>%  count(STR_ID) %>% filter(n==65) %>% dim #303146
trgt_qc_rmlist %>% filter(!(is.na(MC))) %>% count(TRID) %>% filter(n==130) %>% dim() #315767


eh_qc_rmlist %>% filter(FILTER == "PASS") %>% count(STR_ID) %>% filter(n==65) -> eh_qc_rmlist_count
trgt_qc_rmlist %>% filter(!(is.na(MC))) %>% count(TRID) %>% filter(n==130) -> trgt_qc_rmlist_count

head(eh_qc_rmlist_count)
head(trgt_qc_rmlist_count)

VennDiagram::venn.diagram(x = list(long = trgt_qc_rmlist_count$TRID,
                                   short = eh_qc_rmlist_count$STR_ID), filename = NULL,
                          #category.names = c("STR DB", "Long-read", "Short-read"),
                          #col = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원의 색상을 반투명하게 설정
                          #fill = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원 내부 색상도 반투명하게 설정
                          alpha = 0.3, # 원의 반투명도를 설정
                          cex = 1.9, # 카테고리 레이블 크기
                          cat.cex = 1.9, # 집합 이름의 글자 크기
                          cat.fontface = "bold", # 글자 폰트 스타일 설정
                          cat.fontfamily = "Arial") -> p
cowplot::plot_grid(p)

head(eh_patho_qc)
head(ref)
eh_patho_qc %>% filter(FILTER == "PASS") %>% filter(STR_ID %in% ref$ID) %>% count(STR_ID) %>% filter(n==65) %>% dim #73
eh_patho_qc %>% filter(FILTER == "PASS") %>% filter(TYPE == "SPANNING/SPANNING") %>%  count(STR_ID) %>% filter(n==65) %>% dim #66
trgt_patho_qc %>% filter(!(is.na(MC))) %>% count(TRID) %>% filter(n==130) %>% dim() #72  2
trgt_patho_qc %>% filter((is.na(MC)))
trgt_patho_qc %>% count(ID)
#dim(ref)
eh_qc_rmlist %>% filter(FILTER == "PASS") %>% count(STR_ID) %>% filter(n==65) -> eh_qc_rmlist_count
trgt_qc_rmlist %>% filter(!(is.na(MC))) %>% count(TRID) %>% filter(n==130) -> trgt_qc_rmlist_count

eh_qc_rmlist %>% filter(FILTER == "PASS") %>% count(STR_ID) -> eh_qc_rmlist_count
trgt_qc_rmlist %>% filter(!(is.na(MC))) %>% count(TRID) -> trgt_qc_rmlist_count

#eh pathon genimcir region : + 10
head(ref)
head(ru)
ru %>% filter(!str_detect(ID,'chr')) %>% filter(!(ID %in% ref$ID))
ru %>% filter(!str_detect(ID,'chr')) %>% filter(!(ID2 %in% ref$ID2))
# patho 73개
head(eh_qc_rmlist)
eh_qc_rmlist %>% select(STR_ID) %>% unique()%>% dim() #321253
trgt_qc_rmlist %>% select(TRID) %>% unique()%>% dim() #321253


eh_patho_qc %>% filter(STR_ID %in% ref$ID) %>% filter(FILTER == "PASS") %>% count(STR_ID)  -> eh_patho_qc_count
trgt_patho_qc %>% filter(!(is.na(MC))) %>% count(TRID) -> trgt_patho_qc_count
trgt_patho_qc %>% filter(!(is.na(MC))) %>% count(MC == ".")

trgt_qc_rmlist %>% filter(!(is.na(MC))) %>% count(ID) -> trgt_qc_rmlist_sampleCount
eh_qc_rmlist %>% filter(FILTER == "PASS") %>% count(ID) -> eh_qc_rmlist_sampleCount
eh_qc_rmlist %>% filter(FILTER == "PASS") %>% filter(type == "SPANNING/SPANNING") %>% count(ID) -> eh_qc_rmlist_onlyspan_sampleCount
head(eh_qc_rmlist)
eh_patho_qc %>% filter(STR_ID %in% ref$ID) %>% filter(FILTER == "PASS") %>% count(ID)  -> eh_patho_qc_sampleCount
eh_patho_qc %>% filter(STR_ID %in% ref$ID) %>% filter(FILTER == "PASS") %>% filter(TYPE == "SPANNING/SPANNING") %>% count(ID)  -> eh_patho_qc_onlyspan_sampleCount
trgt_patho_qc %>% filter(!(is.na(MC))) %>% count(ID) -> trgt_patho_qc_sampleCount

#trgt_qc_rmlist_sampleCount %>% merge(trgt_patho_qc_sampleCount,by="ID") %>% mutate(prop = (n.x+n.y)/2/321253*100) -> trgt_ID_accu
#eh_qc_rmlist_sampleCount %>% merge(eh_patho_qc_sampleCount,by="ID") %>% mutate(prop = (n.x+n.y)/321253*100) -> eh_ID_accu

#eh_qc_rmlist_onlyspan_sampleCount %>% merge(eh_patho_qc_onlyspan_sampleCount,by="ID") %>% mutate(prop = (n.x+n.y)/321253*100) -> eh_ID_accu_v2

trgt_qc_rmlist_sampleCount %>% merge(trgt_patho_qc_sampleCount,by="ID") %>% mutate(prop = (n.x+n.y)/2/321326*100) -> trgt_ID_accu
eh_qc_rmlist_sampleCount %>% merge(eh_patho_qc_sampleCount,by="ID") %>% mutate(prop = (n.x+n.y)/321326*100) -> eh_ID_accu

eh_qc_rmlist_onlyspan_sampleCount %>% merge(eh_patho_qc_onlyspan_sampleCount,by="ID") %>% mutate(prop = (n.x+n.y)/321326*100) -> eh_ID_accu_v2

trgt_ID_accu <- read.table("~/Desktop/KU/@research/STR/figure/figure1.trgt_pass_prop.v2.txt",header = T)
trgt_ID_accu %>% mutate(prop =  (n.x+n.y)/321326/2*100) -> trgt_ID_accu

eh_ID_accu <- read.table("~/Desktop/KU/@research/STR/figure/figure1.eh_pass_prop.v2.txt",header = T)
eh_ID_accu_v2 <- read.table("~/Desktop/KU/@research/STR/figure/figure1.eh_pass_prop_onlySpan.v2.txt",header = T)
eh_ID_accu %>% mutate(prop =  (n.x+n.y)/321326*100) -> eh_ID_accu
eh_ID_accu_v2 %>% mutate(prop =  (n.x+n.y)/321326*100) -> eh_ID_accu_v2

#write.table(trgt_ID_accu,"~/Desktop/KU/@research/STR/figure/figure1.trgt_pass_prop.v2.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(eh_ID_accu,"~/Desktop/KU/@research/STR/figure/figure1.eh_pass_prop.v2.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(eh_ID_accu_v2,"~/Desktop/KU/@research/STR/figure/figure1.eh_pass_prop_onlySpan.v2.txt",col.names = T,row.names = F,quote = F,sep = "\t")

trgt %>% filter(MC == '0') %>% select(TRID) %>% unique() %>% dim

trgt_ID_accu$type <- "long"
eh_ID_accu$type <- "short"
eh_ID_accu_v2$type <- "short"

trgt_ID_accu %>% rbind(eh_ID_accu) %>%
  ggplot(aes(x=type,y=prop)) + 
  geom_violin() +
  labs(title="Sample Accuracy: Long (PASS) vs Short (PASS)")

trgt_ID_accu %>% rbind(eh_ID_accu_v2) %>%
  ggplot(aes(x=type,y=prop)) + 
  geom_violin() + 
  labs(title="Sample Accuracy: Long (PASS) vs Short (PASS + onlySpanning)")


### suplple
head(eh_qc_rmlist)
head(eh_patho_qc)
eh_patho_qc %>% rename(type = TYPE) ->eh_patho_qc
dim(eh_patho_qc)
head(eh_qc_rmlist)
eh_qc_rmlist %>% rbind(eh_patho_qc %>% filter(STR_ID %in% ref$ID)) %>% count(ID,FILTER,type) -> eh_qc_type_count
#head(eh_patho_qc)
#write.table(eh_qc_type_count,"~/Desktop/KU/@research/STR/figure/sup.figure/eh_qc_alleletype_count.bysample.txt",col.names = T,row.names = F,quote = F,sep = "\t")
trgt_qc_rmlist %>% rbind(trgt_patho_qc) %>% group_by(ID,FILTER,type) %>% count(TRID) -> trgt_qc_type
head(trgt_qc_rmlist)
head(eh_qc_type_count)
table(trgt_patho_qc$TRID) %>% dim()

### venn
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header = T)

head(final_ref)
head(ru)

eh_qc_rmlist %>% rbind(eh_patho_qc %>% filter(STR_ID %in% ref$ID)) %>% filter(FILTER == "PASS",type == "SPANNING/SPANNING") %>% count(STR_ID) -> eh_pass_span_strcount
eh_qc_rmlist %>% rbind(eh_patho_qc %>% filter(STR_ID %in% ref$ID)) %>% filter(FILTER == "PASS") %>% count(STR_ID) -> eh_pass_strcount
trgt_qc_rmlist %>% rbind(trgt_patho_qc) %>% filter(!(is.na(MC))) %>% count(TRID) -> trgt_pass_strcount
setwd("~/Desktop/KU/@research/STR/figure/")
head(trgt_pass_strcount)
trgt_pass_strcount %>% filter(n==130) -> trgt_pass_strcount_allsample
eh_pass_span_strcount %>% filter(n==65) -> eh_pass_span_strcount_allsample
eh_pass_strcount %>% filter(n==65) -> eh_pass_strcount_allsample
head(trgt_pass_strcount_allsample)
trgt_pass_strcount_allsample %>% select(TRID) %>% rename(TRGT = TRID) %>% mutate(STR_ID = TRGT) -> trgt_ID
eh_pass_strcount_allsample %>% select(STR_ID) %>% rename(EH = STR_ID) %>% mutate(STR_ID = EH) -> eh_ID

str(final_ref$ID)
venndata <- list(
  db = final_ref$ID,
  long = trgt_pass_strcount_allsample$TRID,
  short=eh_pass_strcount_allsample$STR_ID
    #short=eh_pass_span_strcount_allsample$STR_ID
)
final_ref %>% select(ID) %>% rename(STR_ID = ID) %>% left_join(trgt_ID) %>% left_join(eh_ID) %>% rename(STR_DB = STR_ID)%>%
  write.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",row.names = F,col.names = T,quote = F,sep = "\t")

#ru %>% select(ID) %>% rename(STR_ID = ID) %>% left_join(trgt_ID) %>% left_join(eh_ID) %>%
#rename(STR_DB = STR_ID) %>%
#write.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_STR_ID_forVenn.txt",row.names = F,col.names = T,quote = F,sep = "\t")


VennDiagram::venn.diagram(x = venndata, filename = NULL,
                          category.names = c("STR DB", "Long-read", "Short-read"),
                          col = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원의 색상을 반투명하게 설정
                          fill = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원 내부 색상도 반투명하게 설정
                          alpha = 0.3, # 원의 반투명도를 설정
                          cex = 1.9, # 카테고리 레이블 크기
                          cat.cex = 1.9, # 집합 이름의 글자 크기
                          cat.fontface = "bold", # 글자 폰트 스타일 설정
                          cat.fontfamily = "Arial") -> p
cowplot::plot_grid(p)


##### 

venn_rawdata %>% filter(is.na(EH)) %>% filter(!is.na(TRGT))
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR
head(common_STR);dim(common_STR) #311542
final_ref %>% filter(ID %in% common_STR$STR_DB) %>% select(MOTIFS) -> common_motif

common_motif %>% mutate(RU.length = ifelse(str_detect(MOTIFS,","),"CRU",str_length(MOTIFS))) %>%
  count(RU.length) -> length_common_motif

#write.table(length_common_motif,"figure1.length_common_motif.count.txt",col.names = T,row.names = F,quote = F,sep = "\t")
ru.scale = c(seq(2,19),"20+")
ru.scale1 = c(seq(2,6))
ru.scale2 = c(seq(7,19),"20+")

seq(2,5)
#table(length_common_motif$length)
head(length_common_motif)
#length_common_motif %>% is.na()
length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(facet = ifelse(RU.length %in% ru.scale1,1,2)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Length of Repeat Units") + 
  facet_wrap(~facet,scales = "free") + 
  theme(strip.text = element_blank())


length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>%
  filter(RU.length %in% ru.scale1) %>% 
  ggplot(aes(x=factor(RU.length,levels=ru.scale1),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  theme(axis.title.x = element_blank()) + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs") -> p1

#x="Length of Repeat Units"
p1

length_common_motif %>% #filter(is.na(length))
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>%
  filter(RU.length %in% ru.scale2) %>% 
  ggplot(aes(x=factor(RU.length,levels=ru.scale2),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) -> p2

combined_plot <- plot_grid(p1)
final_plot <- add_sub(combined_plot, 
                      "Common X-axis Title", 
                      vjust = 0.5, 
                      size = 14, 
                      fontface = "bold",    # 텍스트 굵게 설정
                      hjust = 0.5)      
cowplot::plot_grid(final_plot)


# 첫 번째 플롯 생성
plot1 <- length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale, RU.length, "20+")) %>% 
  mutate(facet = ifelse(RU.length %in% ru.scale1, 1, 2)) %>%
  filter(facet == 1) %>%
  ggplot(aes(x = factor(RU.length, levels = ru.scale), y = n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) +
  theme(axis.title.x= element_blank())
plot2 <- length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale, RU.length, "20+")) %>% 
  mutate(facet = ifelse(RU.length %in% ru.scale1, 1, 2)) %>%
  filter(facet == 2) %>%
  ggplot(aes(x = factor(RU.length, levels = ru.scale), y = n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) +
  theme(axis.title.x= element_blank(),
        axis.title.y= element_blank())

combined_plot <- plot_grid(plot1, plot2, nrow = 1, rel_widths = c(1, 2))

# 공통 xlab, ylab 추가 및 theme 적용
final_plot <- ggdraw() +
  draw_plot(combined_plot, 0, 0, 1, 1) +
  draw_label("Length of Repeat Units", x = 0.5, y = 0, vjust = 0, size = 16) +
  draw_label("# of STRs", x = 0, y = 0.5, angle = 90, vjust = 1, size = 16) + 
  theme_step1()
plot_grid(final_plot)




length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(facet = ifelse(RU.length %in% ru.scale1,1,2)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Length of Repeat Units") + 
  facet_row(~ facet, scales = "free", space = "free") + 
  theme(strip.text = element_blank()) -> ep2.1





############ STR processing
load(file="~/Desktop/KU/@research/STR/02.compare/STR_prep/eh.ganstr.v13.raw.merge.RDATA")
load(file="~/Desktop/KU/@research/STR/02.compare/STR_prep/trgt.ganstr.v13.raw.merge.RData")

eh_patho <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_prep/eh.pathogenic_analysisDB.raw.merge.txt",header = T)
eh_patho %>% rename(type = TYPE) ->eh_patho
trgt_patho <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_prep/trgt.pathogenic_analysisDB.raw.merge.txt",header = T)
ref <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/db/str-analysis/str_analysis/variant_catalogs/catalog.GRCh38.with_adjacent_repeats.TRGT.bed",col_names = F)
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)

ref %>%  mutate(MOTIFS = str_split_fixed(X4,";",3)[,2]) %>% 
  mutate(MOTIFS = str_split_fixed(MOTIFS,"=",2)[,2]) %>% #head()
  mutate(ID = str_split_fixed(X4,";",3)[,1]) %>% 
  mutate(ID = str_split_fixed(ID,"=",2)[,2]) %>% #head()
  select(X1,X2,X3,MOTIFS,ID) -> ref


#final_ref %>% filter(str_detect(ID, "ARX"))
head(eh_patho)

head(eh)
head(trgt)
head(final_ref)
#head(rmlist_final)

head(trgt_patho)
trgt_patho %>% filter(str_detect(MOTIFS,",")) %>% select(TRID) %>% unique() -> complex_str_list
#eh %>% filter((str_detect(STR_ID,"chr"))) %>% dim()

eh %>% filter((str_detect(STR_ID,"chr"))) %>% filter(STR_ID %in% final_ref$ID) %>% filter(FILTER == "PASS") -> eh_simple_qc
pattern <- str_c(complex_str_list$TRID, collapse = "|")
head(eh_patho)
eh_patho %>% filter(str_detect(STR_ID, pattern)) %>% select(STR_ID,RU) %>% unique() %>% mutate(TRID = str_split_fixed(STR_ID,"_",2)[,1]) -> complex_str_list_forEH
#write.table(complex_str_list_forEH,"~/Desktop/KU/@research/STR/figure/extra_info/complex_STR.v2.txt",row.names = F,col.names = T,quote = F,sep = "\t")

head(complex_str_list_forEH)

eh_patho %>% filter(!(STR_ID %in% complex_str_list_forEH$STR_ID)) -> eh_patho_simple_qc
eh_patho %>% filter(STR_ID %in% complex_str_list_forEH$STR_ID) -> eh_patho_complex_qc
head(venn_rawdata)
colnames(eh_patho_simple_qc)[1:2] <- c("CHROM","POS")
head(eh_patho_complex_qc)
colnames(eh_patho_complex_qc)[1:2] <- c("CHROM","POS")
eh_complex_pass <- eh_patho_complex_qc

eh_simple_qc %>% rbind(eh_patho_simple_qc) %>% filter(STR_ID %in% venn_rawdata$EH) -> eh_simple_pass
eh_simple_pass %>% filter(STR_ID %in% venn_rawdata$TRGT) -> eh_simple_pass_intersect
eh_simple_pass %>% filter(!(STR_ID %in% venn_rawdata$TRGT)) -> eh_simple_pass_onlyEH


#/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2
#save(eh_simple_pass_intersect, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_simple_pass_intersect.RData")
#save(eh_simple_pass_onlyEH, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_simple_pass_onlyEH.RData")
#save(eh_complex_pass, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_complex_pass.RData")

#rm(eh_simple_qc)
head(final_ref)
final_ref %>% filter(!(str_detect(ID,"chr")))

trgt %>% filter((str_detect(TRID,"chr"))) %>% filter(TRID %in% final_ref$ID) %>% filter(!(is.na(MC))) -> trgt_simple_qc
trgt_patho %>% filter(TRID %in% complex_str_list_forEH$TRID) %>% filter(!(is.na(MC))) -> trgt_patho_complex_qc
trgt_patho %>% filter(!(TRID %in% complex_str_list_forEH$TRID)) %>% filter(!(is.na(MC))) -> trgt_patho_simple_qc
head(trgt_patho_simple_qc)
table(trgt_patho_complex_qc$TRID) %>% dim()
#head(trgt_patho)
head(complex_str_list_forEH)
trgt_complex_pass <-trgt_patho_complex_qc


trgt_simple_qc %>% rbind(trgt_patho_simple_qc) %>% filter(TRID %in% venn_rawdata$TRGT) -> trgt_simple_pass
trgt_simple_pass %>% filter(TRID %in% venn_rawdata$EH) -> trgt_simple_pass_intersect
trgt_simple_pass %>% filter(!(TRID %in% venn_rawdata$EH)) -> trgt_simple_pass_onlytrgt

save(trgt_simple_pass_intersect, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_simple_pass_intersect.RData")
save(trgt_simple_pass_onlytrgt, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_simple_pass_onlytrgt.RData")
save(trgt_complex_pass, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_complex_pass.RData")


head(trgt_simple_pass_intersect)



###### STR length : only EH vs only TRGT
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR
head(venn_rawdata)
venn_rawdata %>% filter(is.na(EH)) %>% select(-EH) %>% na.omit() %>% dim()

final_ref %>% filter(ID %in% venn_rawdata$TRGT) %>% filter(!(ID %in% venn_rawdata$EH)) %>% dim() #4297
final_ref %>% filter(ID %in% venn_rawdata$TRGT) %>% filter(!(ID %in% venn_rawdata$EH)) -> onlyTRGT
final_ref %>% filter(ID %in% venn_rawdata$EH) %>% filter(!(ID %in% venn_rawdata$TRGT)) -> onlyEH


ru.scale = c(seq(2,19),"20+")
ru.scale1 = c(seq(2,6))
ru.scale2 = c(seq(7,19),"20+")


onlyTRGT %>% head()
library(ggforce)
onlyEH$type <- "SRS"
onlyTRGT$type <- "LRS"

onlyEH %>% rbind(onlyTRGT) %>% group_by(type) %>%
  mutate(RU.length = str_length(MOTIFS)) %>% count(RU.length) %>% 
  write.table("~/Desktop/KU/@research/STR/figure/figure1.length_unique_motif.count.LvsS.txt",row.names = F,col.names = T,quote = F,sep = "\t")

onlyEH %>% rbind(onlyTRGT) %>% group_by(type) %>%
  mutate(RU.length = str_length(MOTIFS)) %>% count(RU.length) %>% #head()
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(facet = ifelse(RU.length %in% ru.scale1,1,2)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n,fill=type)) + 
  geom_bar(stat = 'identity',position = 'dodge') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Length of Repeat Units") + 
  facet_row(~ facet, scales = "free", space = "free") + 
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme(strip.text = element_blank())



  
############################################################################## STR unit count
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)

venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR
head(final_ref)
head(common_STR)
final_ref %>% filter(ID %in% venn_rawdata$TRGT) %>% filter(!(ID %in% venn_rawdata$EH)) %>% dim() #4297
final_ref %>% filter(ID %in% venn_rawdata$TRGT) %>% filter(!(ID %in% venn_rawdata$EH)) %>% select(MOTIFS) -> onlyTRGT
final_ref %>% filter(ID %in% venn_rawdata$EH) %>% filter(!(ID %in% venn_rawdata$TRGT)) %>% select(MOTIFS) -> onlyEH

onlyEH$type <- "SRS"
onlyTRGT$type <- "LRS"



final_ref %>% filter(ID %in% common_STR$STR_DB) %>% select(MOTIFS) -> common_motif

ru.scale = c(seq(1,19),"20+")
ru.scale1 = c(seq(1,2))
ru.scale2 = c(seq(3,19),"20+")

onlyTRGT %>% rbind(onlyEH) %>% group_by(type) %>% count(MOTIFS) %>% arrange(-n) %>% #mutate(Rank = rank(-n)) %>% #head()
  mutate(new_n = ifelse(n %in% c(1:19),n,"20+")) %>% #count(new_n)
  count(new_n) %>% #filter(!(new_n %in% ru.scale))
  mutate(facet = ifelse(new_n %in% ru.scale1,1,2)) %>%
  ggplot(aes(x=factor(new_n,levels=ru.scale),y=n,fill=type)) + 
  geom_bar(stat = 'identity',position = 'dodge') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(x="Frequency of Repeat Units",y="# of STRs") +
  facet_row(~ facet, scales = "free", space = "free") + 
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme(strip.text = element_blank())

common_motif  %>% count(MOTIFS) %>% arrange(-n) %>% #mutate(Rank = rank(-n)) %>% #head()
  mutate(new_n = ifelse(n %in% c(1:19),n,"20+")) %>% #count(new_n)
  count(new_n) %>% #filter(!(new_n %in% ru.scale))
  mutate(facet = ifelse(new_n %in% ru.scale1,1,2)) %>%
  ggplot(aes(x=factor(new_n,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity',position = 'dodge') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(x="Frequency of Repeat Units",y="# of STRs") +
  facet_row(~ facet, scales = "free", space = "free") + 
  theme(strip.text = element_blank()) ->ep2.2
#ep2.2

#ep2.1,ep2.2
ggarrange(ep2.1,ep2.2,ncol = 1,labels = c('A',"B"), font.label = list(size = 28), label.y = 1.01) -> ep2
ep2

####### no calling length (motif lengg * repeat count)


load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_simple_pass_onlytrgt.RData")
load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_simple_pass_onlyEH.RData")

head(trgt_simple_pass_onlytrgt)
head(eh_simple_pass_onlyEH)

eh_simple_pass_onlyEH %>% select(-type,-ADSP,-ADFL,-ADIR) %>% #filter(FILTER == "PASS") %>%
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>%
  mutate(ALT = ifelse(ALT == ".",0,ALT)) %>% mutate(GT = str_split_fixed(GT,"/",2),ALT= str_split_fixed(ALT,",",2)) -> eh_simple_pass_onlyEH_gt

eh_simple_pass_onlyEH_gt %>% mutate(STR1 = ifelse(GT[,1] == "1",ALT[,1],ifelse(GT[,1] == "2",ALT[,2],REF))) %>%
  mutate(STR2 = ifelse(GT[,2] == "1",ALT[,1],ifelse(GT[,2] == "2",ALT[,2],REF))) %>% select(-ALT,-GT,-REF) %>% mutate(STR1 = as.integer(STR1),STR2 = as.integer(STR2)) -> eh_simple_pass_onlyEH_gt

head(eh_simple_pass_onlyEH)

eh_simple_pass_onlyEH %>% select(ID,STR_ID,type) %>% mutate(type1 = str_split_fixed(type,"/",2)[,1],type2 = str_split_fixed(type,"/",2)[,2]) %>% 
  select(-type) -> eh_simple_pass_onlyEH_type

head(eh_simple_pass_onlyEH_gt)
head(eh_simple_pass_onlyEH_type)

eh_simple_pass_onlyEH_gt %>% select(-FILTER) %>% left_join(eh_simple_pass_onlyEH_type) -> eh_simple_pass_onlyEH_processing

eh_simple_pass_onlyEH_processing %>% mutate(EH_STR1 = ifelse(STR1 < STR2, STR1,STR2),EH_STR2 = ifelse(STR1 < STR2, STR2,STR1)) %>%
  mutate(EH_type1 = ifelse(STR1 < STR2, type1,type2),EH_type2 = ifelse(STR1 < STR2, type2,type1)) %>% 
  select(-STR1,-STR2,-type1,-type2)-> eh_simple_pass_onlyEH_processing

head(eh_simple_pass_onlyEH_processing)

eh_simple_pass_onlyEH_processing %>% mutate(RU.length = str_length(RU)) %>% #head()
  mutate(EH_STR1_length = EH_STR1*RU.length,EH_STR2_length = EH_STR2*RU.length) -> eh_simple_pass_onlyEH_processing

head(a)
eh_simple_pass_onlyEH_processing %>%
  select(STR_ID,ID,EH_type1,EH_type2) %>% 
  pivot_longer(cols = c("EH_type1","EH_type2")) -> a

eh_simple_pass_onlyEH_processing %>%
  select(STR_ID,ID,EH_STR1_length,EH_STR2_length) %>%
  pivot_longer(cols = c("EH_STR1_length","EH_STR2_length")) -> b

head(a)  
head(b)
a %>% rename(type = value,allele = name) %>% mutate(allele = ifelse(allele == "EH_type1","allele1","allele2")) -> a
b %>% rename(STR_length = value,allele = name) %>% mutate(allele = ifelse(allele == "EH_STR1_length","allele1","allele2")) ->b

a %>% left_join(b) %>%
  ggplot(aes(x=type,y=STR_length,fill=type)) + 
  geom_violin()

head(b)
pivot_longer(cols = c("EH_type1","EH_type2")) + 
  ggplot(aes(x=EH_STR1_length,y=EH_STR2_length)) + 
  geom_point()
  

##
head(eh)
table(eh$ID) %>% dim
head(trgt)
#eh %>% count(FILTER)
head(venn_rawdata)
venn_rawdata %>% filter(is.na(TRGT) | is.na(EH)) %>% filter(!(is.na(TRGT) & is.na(EH))) -> not_intersect
head(not_intersect)

eh %>% filter(STR_ID %in% not_intersect$STR_DB) %>% filter(ID != "NIH20N2594890") -> eh_notintersect
trgt %>% filter(TRID %in% not_intersect$STR_DB) %>% filter(ID != "NIH23J3904558") -> trgt_notintersect

#eh_notintersect %>% count(FILTER)

head(trgt_notintersect)
trgt_patho %>%  filter(TRID %in% not_intersect$STR_DB) %>% filter(ID != "NIH23J3904558") -> trgt_patho_notintersect
trgt_notintersect %>% rbind(trgt_patho_notintersect) -> trgt_notintersect

head(trgt_notintersect)


trgt_notintersect
head(eh_notintersect)
eh_notintersect %>% select(-type,-ADSP,-ADFL,-ADIR) %>% #filter(FILTER == "PASS") %>%
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>% #head()
  mutate(ALT = ifelse(ALT == ".",0,ALT)) %>% #head()
  mutate(GT = str_split_fixed(GT,"/",2),ALT= str_split_fixed(ALT,",",2)) -> eh_notintersect_gt
head(eh_notintersect_gt)

eh_notintersect_gt %>% mutate(STR1 = ifelse(GT[,1] == "1",ALT[,1],ifelse(GT[,1] == "2",ALT[,2],REF))) %>%
  mutate(STR2 = ifelse(GT[,2] == "1",ALT[,1],ifelse(GT[,2] == "2",ALT[,2],REF))) %>% select(-ALT,-GT,-REF) %>% mutate(STR1 = as.integer(STR1),STR2 = as.integer(STR2)) -> eh_notintersect_gt

head(eh_notintersect_gt)

eh_notintersect %>% select(ID,STR_ID,type) %>% mutate(type1 = str_split_fixed(type,"/",2)[,1],type2 = str_split_fixed(type,"/",2)[,2]) %>% #head()
  select(-type) -> eh_notintersect_type

head(eh_notintersect_gt)
head(eh_notintersect_type)

eh_notintersect_gt %>%  left_join(eh_notintersect_type) -> eh_notintersect_processing

eh_notintersect_processing %>% mutate(EH_STR1 = ifelse(STR1 < STR2, STR1,STR2),EH_STR2 = ifelse(STR1 < STR2, STR2,STR1)) %>%
  mutate(EH_type1 = ifelse(STR1 < STR2, type1,type2),EH_type2 = ifelse(STR1 < STR2, type2,type1)) %>% 
  select(-STR1,-STR2,-type1,-type2)-> eh_notintersect_processing

eh_notintersect_processing %>% select(-CHROM,-POS) -> eh_notintersect_processing

head(eh_notintersect_processing)
head(trgt_notintersect)
head(trgt_notintersect)
#trgt_notintersect %>% filter(is.na(AP))

trgt_notintersect %>% filter(!is.na(AP)) %>% select(-CHROM,-POS,-STRUC) -> trgt_notintersect


head(trgt_notintersect)




#write.table(eh_notintersect_processing,"~/Desktop/KU/@research/STR/figure/figure1.eh_notintersect_processing.txt",col.names = T,row.names = F,quote = F,sep="\t")
#write.table(trgt_notintersect,"~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect.txt",col.names = T,row.names = F,quote = F,sep="\t")


##
eh_notintersect_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure1.eh_notintersect_processing.txt")
trgt_notintersect <- read_table("~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect.txt")

head(trgt_notintersect)

head(eh_notintersect_processing)
head(trgt_notintersect)

trgt_notintersect %>% select(ID,TRID,MOTIFS,Allele,MC) %>% #head()
  pivot_wider(names_from = Allele,values_from = MC) %>% rename("MC1" = allele_1,"MC2" = allele_2) -> trgt_notintersect_gt

trgt_notintersect %>% select(ID,TRID,MOTIFS,Allele,AP) %>% #head()
  pivot_wider(names_from = Allele,values_from = AP) %>% rename("AP1" = allele_1,"AP2" = allele_2) -> trgt_notintersect_AP

trgt_notintersect %>% select(ID,TRID,MOTIFS,Allele,AM) %>% #head()
  pivot_wider(names_from = Allele,values_from = AM) %>% rename("AM1" = allele_1,"AM2" = allele_2) -> trgt_notintersect_AM


head(trgt_notintersect_gt)
head(trgt_notintersect_AM)

trgt_notintersect_gt %>% left_join(trgt_notintersect_AM) %>% left_join(trgt_notintersect_AP) -> trgt_notintersect_merge
head(trgt_notintersect_merge)
head(trgt_notintersect_merge)
trgt_notintersect_merge %>% filter(MC1 == ".")
trgt_notintersect_merge %>% filter(MC2 == ".")

trgt_notintersect_merge %>% mutate(across(AM1:AP2, ~na_if(., "."))) ->  trgt_notintersect_merge

head(eh_notintersect_processing)
head(trgt_notintersect_merge)

trgt_notintersect_merge %>% mutate(TRGT_STR1 = ifelse(MC1 < MC2, MC1,MC2),TRGT_STR2 = ifelse(MC1 < MC2, MC2,MC1)) %>%
  mutate(TRGT_AM1 = ifelse(MC1 < MC2, AM1,AM2),TRGT_AM2 = ifelse(MC1 < MC2, AM2,AM1)) %>%
  mutate(TRGT_AP1 = ifelse(MC1 < MC2, AP1,AP2),TRGT_AP2 = ifelse(MC1 < MC2, AP2,AP1)) %>% 
  select(-MC1,-MC2,-AM1,-AM2,-AP1,-AP2)  -> trgt_notintersect_merge_processing

#write.table(trgt_notintersect_merge_processing,"~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect_merge_processing.txt",col.names = T,row.names = F,quote = F,sep="\t")


head(eh_notintersect_processing)
head(trgt_notintersect_merge_processing)

### call rate

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
#head(final_ref)
head(venn_rawdata)

venn_rawdata %>% filter(is.na(TRGT),is.na(EH)) -> onlyDB
head(onlyDB)
venn_rawdata %>% filter(is.na(EH),!is.na(TRGT))  -> onlyTRGT
head(onlyTRGT)
venn_rawdata %>% filter(!is.na(EH),is.na(TRGT))  -> onlyEH
head(onlyEH)


eh_notintersect_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure1.eh_notintersect_processing.txt")
trgt_notintersect <- read_table("~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect.txt")
head(eh_notintersect_processing)
head(trgt_notintersect)
eh_notintersect_processing %>% filter(FILTER == "PASS") %>% count(STR_ID) %>%  mutate(EH_pass=n/65) %>% select(-n) -> eh_notintersect_processing_callrate
trgt_notintersect %>% select(ID,TRID,Allele,MC) %>% count(TRID) %>% mutate(TRGT_pass=n/130) %>% rename(STR_ID = TRID) %>% select(-n) -> trgt_notintersect_callrate


head(final_ref)
final_ref %>% mutate(RU.length = str_length(MOTIFS)) %>% select(MOTIFS,ID,RU.length) %>% rename(STR_ID = ID)-> final_ref_rulengt
head(final_ref_rulengt)
#final_ref %>% filter(!str_detect(ID,"chr"))

head(eh_notintersect_processing_callrate)
head(trgt_notintersect_callrate)
#head(final_ref)
eh_notintersect_processing_callrate %>% left_join(trgt_notintersect_callrate) %>% mutate(g = ifelse(STR_ID %in% onlyTRGT$STR_DB,"TRGT",ifelse(STR_ID %in% onlyEH$STR_DB,"EH","none"))) %>% #head()
  left_join(final_ref_rulengt) -> merge_call_rate

head(merge_call_rate)
merge_call_rate %>% count(g)
head(final_ref_rulengt)
##

eh_notintersect_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure1.eh_notintersect_processing.txt")
trgt_notintersect_merge_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect_merge_processing.txt")

head(merge_call_rate)
head(eh_notintersect_processing)
head(trgt_notintersect_merge_processing)
merge_call_rate %>% filter(STR_ID == "chr11_96563250_96563265")

trgt_notintersect_merge_processing %>% filter(TRID == "chr11_96563250_96563265")
head(final_ref)
venn_rawdata %>% filter(STR_DB == "chr10_100321663_100321755")
eh_notintersect_processing %>% select(STR_ID,ID,EH_STR1,EH_STR2) %>% filter(STR_ID %in% onlyEH$EH) %>% head()
head(onlyTRGT)

venn_rawdata %>% filter(STR_DB == "chr10_118248242_118248314")
#filter(TRID == "chr10_118248242_118248314")
trgt_notintersect_merge_processing %>% filter(TRID %in% onlyTRGT$TRGT) %>% #filter(is.na(TRGT_AP1))
  pivot_longer(TRGT_STR1:TRGT_AP2) %>% 
  mutate(name = ifelse(str_detect(name,"STR"),"MC",ifelse(str_detect(name,"AM"),"AM","AP"))) %>%
  group_by(TRID,name) %>% #head()
  summarise(value = mean(value)) %>% 
  pivot_wider(names_from = name,values_from = value) -> tmp
  
head(tmp)
head(final_ref_rulengt)
tmp %>% rename(STR_ID = TRID) %>% left_join(final_ref_rulengt) %>% left_join(merge_call_rate %>% select(-TRGT_pass,MOTIFS,RU.length)) %>%  
  ungroup() %>% mutate(EH_pass = ifelse(is.na(EH_pass),0,EH_pass)) %>% mutate(STR_length = MC * RU.length) %>% select(-MC,-g) -> tmp

tmp %>% filter(is.na(AM))
tmp %>% filter(is.na(AP))
tmp %>% filter(is.na(AP))
venn_rawdata %>% filter(STR_DB == "chr10_118248242_118248314")
trgt_notintersect_merge_processing %>% filter(TRID == "chr10_118248242_118248314")

head(tmp)
tmp %>% mutate(AM = ifelse(is.na(AM),-1,AM)) %>% filter(is.na(AP))
tmp %>% ggplot(aes(x=AP,y=AM)) + 
  geom_point() + 
  theme_step1()



tmp %>% ggplot(aes(x=STR_length,y=EH_pass,fill= AP)) + 
  geom_point()

tmp %>% ggplot(aes(x=STR_length,y=AM,fill= AP)) + 
  geom_point()
head(trgt_notintersect_merge_processing)

trgt_notintersect_merge_processing %>% filter(is.na(TRGT_AP1)) %>%  count(TRGT_STR1,TRGT_AP1)
trgt_notintersect_merge_processing %>% filter(is.na(TRGT_AP2)) %>%  count(TRGT_STR2,TRGT_AP2)

trgt_notintersect_merge_processing %>% filter(TRGT_STR1 == 0) %>%  count(TRGT_STR1,TRGT_AM1)
trgt_notintersect_merge_processing %>% filter(is.na(TRGT_STR2)) %>%  count(TRGT_STR2,TRGT_AM2)

trgt_notintersect_merge_processing %>% ggplot(aes(x=TRGT_AM1,y=TRGT_AM2)) + geom_point()
  

trgt_notintersect_merge_processing %>% filter(TRID %in% onlyTRGT$TRGT) %>% select(ID:TRGT_STR2) %>%
  pivot_longer(TRGT_STR1:TRGT_STR2,values_to = "MC") %>% mutate(name = ifelse(str_detect(name,"STR1"),"STR1","STR2")) -> trgt_notintersect_merge_processing_MC

head(trgt_notintersect_merge_processing_MC)
head(trgt_notintersect_merge_processing)
trgt_notintersect_merge_processing %>% filter(TRID %in% onlyTRGT$TRGT) %>% select(ID,TRID,TRGT_AM1,TRGT_AM2) %>%
  pivot_longer(TRGT_AM1:TRGT_AM2,values_to = "AM") %>% mutate(name = ifelse(str_detect(name,"1"),"STR1","STR2")) %>% 
  mutate(AM = ifelse(is.na(AM),0,AM)) -> trgt_notintersect_merge_processing_AM

trgt_notintersect_merge_processing %>% filter(TRID %in% onlyTRGT$TRGT) %>% select(ID,TRID,TRGT_AP1,TRGT_AP2) %>%
  pivot_longer(TRGT_AP1:TRGT_AP2,values_to = "AP") %>% mutate(name = ifelse(str_detect(name,"1"),"STR1","STR2")) %>% 
  mutate(AP = ifelse(is.na(AP),1,AP)) -> trgt_notintersect_merge_processing_AP


head(merge_call_rate)
merge_call_rate %>% count(g)

trgt_notintersect_merge_processing_MC %>% left_join(trgt_notintersect_merge_processing_AM) %>% left_join(trgt_notintersect_merge_processing_AP) %>%
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(STR_length = MC * RU.length) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>% #head()
  rename(STR_ID = TRID) %>% left_join(merge_call_rate %>% select(STR_ID,EH_pass,g)) %>% #head()
  mutate(EH_pass = ifelse(is.na(EH_pass),0,EH_pass))-> tmp2

head(tmp2)

head(tmp2)
tmp2 %>% count(g)
tmp2 %>% filter(is.na(g)) %>% count(EH_pass)
tmp2 %>% mutate(EH_pass = ifelse(is.na(EH_pass,0,EH_pass)))


c(0, seq(0.1, 1, by = 0.1)) 


tmp2 %>% mutate(EH_pass_range = case_when(
  EH_pass == 0 ~ "0",
  EH_pass > 0 & EH_pass < 0.1 ~ "(0,0.1)",
  EH_pass >= 0.1 & EH_pass < 0.2 ~ "[0.1,0.2)",
  EH_pass >= 0.2 & EH_pass < 0.3 ~ "[0.2,0.3)",
  EH_pass >= 0.3 & EH_pass < 0.4 ~ "[0.3,0.4)",
  EH_pass >= 0.4 & EH_pass < 0.5 ~ "[0.4,0.5)",
  EH_pass >= 0.5 & EH_pass < 0.6 ~ "[0.5,0.6)",
  EH_pass >= 0.6 & EH_pass < 0.7 ~ "[0.6,0.7)",
  EH_pass >= 0.7 & EH_pass < 0.8 ~ "[0.7,0.8)",
  EH_pass >= 0.8 & EH_pass < 0.9 ~ "[0.8,0.9)",
  EH_pass >= 0.9 & EH_pass < 1 ~ "[0.9,1)")) -> tmp2

tmp2 %>% filter(EH_pass == 0) %>% count(EH_pass_range)
head(tmp2)

tmp2 %>% count(EH_pass_range)
tmp2 %>% ggplot(aes(x=STR_length,y=GC)) + 
  geom_point() + 
  facet_wrap(~EH_pass_range,ncol = 4)

tmp2 %>% ggplot(aes(x=factor(EH_pass_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                         "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                         "[0.8,0.9)", "[0.9,1)")),y=GC,fill=EH_pass_range)) + 
  geom_violin()

head(tmp2)

tmp2 %>% ggplot(aes(x=STR_length,y=AM)) + 
  geom_point() + 
  facet_wrap(~EH_pass_range,ncol = 4)

tmp2 %>% ggplot(aes(x=factor(EH_pass_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                     "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                     "[0.8,0.9)", "[0.9,1)")),y=AM,fill=EH_pass_range)) + 
  geom_violin()

tmp2 %>% ggplot(aes(x=factor(EH_pass_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                     "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                     "[0.8,0.9)", "[0.9,1)")),y=AP,fill=EH_pass_range)) + 
  geom_violin()

head(tmp2)
tmp2 %>% ggplot(aes(x=as.factor(GC),y=EH_pass,fill=GC)) + 
  geom_violin()


tmp2 %>% select(STR_ID,STR_length,AM,AP,GC,EH_pass) %>%
  group_by(STR_ID) %>% #head()
  summarise(across(STR_length:EH_pass,mean)) %>%
  mutate(EH_pass_range = case_when(
    EH_pass == 0 ~ "0",
    EH_pass > 0 & EH_pass < 0.1 ~ "(0,0.1)",
    EH_pass >= 0.1 & EH_pass < 0.2 ~ "[0.1,0.2)",
    EH_pass >= 0.2 & EH_pass < 0.3 ~ "[0.2,0.3)",
    EH_pass >= 0.3 & EH_pass < 0.4 ~ "[0.3,0.4)",
    EH_pass >= 0.4 & EH_pass < 0.5 ~ "[0.4,0.5)",
    EH_pass >= 0.5 & EH_pass < 0.6 ~ "[0.5,0.6)",
    EH_pass >= 0.6 & EH_pass < 0.7 ~ "[0.6,0.7)",
    EH_pass >= 0.7 & EH_pass < 0.8 ~ "[0.7,0.8)",
    EH_pass >= 0.8 & EH_pass < 0.9 ~ "[0.8,0.9)",
    EH_pass >= 0.9 & EH_pass < 1 ~ "[0.9,1)")) %>%
  pivot_longer(STR_length:GC) %>% #filter(name != "RU.leng")
  ggplot(aes(x=factor(EH_pass_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                              "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                              "[0.8,0.9)", "[0.9,1)")),y=value,fill=factor(EH_pass_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                                                                                   "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                                                                                   "[0.8,0.9)", "[0.9,1)")))) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  facet_wrap(~name,ncol = 2,scales = 'free') + 
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())



tmp2 %>% select(STR_ID,STR_length,AM,AP,GC,EH_pass) %>%
  group_by(STR_ID) %>% #head()
  summarise(across(STR_length:EH_pass,mean)) %>%
  mutate(EH_pass_range = case_when(
    EH_pass == 0 ~ "0",
    EH_pass > 0 & EH_pass < 0.1 ~ "(0,0.1)",
    EH_pass >= 0.1 & EH_pass < 0.2 ~ "[0.1,0.2)",
    EH_pass >= 0.2 & EH_pass < 0.3 ~ "[0.2,0.3)",
    EH_pass >= 0.3 & EH_pass < 0.4 ~ "[0.3,0.4)",
    EH_pass >= 0.4 & EH_pass < 0.5 ~ "[0.4,0.5)",
    EH_pass >= 0.5 & EH_pass < 0.6 ~ "[0.5,0.6)",
    EH_pass >= 0.6 & EH_pass < 0.7 ~ "[0.6,0.7)",
    EH_pass >= 0.7 & EH_pass < 0.8 ~ "[0.7,0.8)",
    EH_pass >= 0.8 & EH_pass < 0.9 ~ "[0.8,0.9)",
    EH_pass >= 0.9 & EH_pass < 1 ~ "[0.9,1)")) %>%
  pivot_longer(STR_length:GC) %>% #filter(name != "RU.leng")
  ggplot(aes(x=EH_pass,y=value)) + 
  geom_point() + 
  facet_wrap(~name,ncol = 2,scales = 'free')


head(eh_notintersect_processing)
head(trgt_notintersect_merge_processing)

eh_notintersect_processing %>% filter(EH_STR1 == 0) %>% filter(FILTER == "PASS")

eh_notintersect_processing %>% filter(STR_ID %in% onlyEH$EH) %>% select(STR_ID,ID,EH_STR1,EH_STR2) %>% #head()
  pivot_longer(EH_STR1:EH_STR2) %>% group_by(STR_ID) %>%
  summarise(STR_mean = mean(value)) %>% mutate(g = "EH") -> a


trgt_notintersect_merge_processing %>% filter(TRID %in% onlyTRGT$TRGT) %>% rename(STR_ID = TRID) %>% select(STR_ID,ID,TRGT_STR1,TRGT_STR2) %>% #head()
  pivot_longer(TRGT_STR1:TRGT_STR2) %>% group_by(STR_ID) %>%
  summarise(STR_mean = mean(value)) %>% mutate(g = "TRGT")-> b

head(a)
head(b)
a %>% rbind(b) %>%dim()#6998
head(final_ref)
head(final_ref_rulengt)
a %>% rbind(b) %>% #left_join(final_ref %>% rename(STR_ID = ID) %>% select(MOTIFS,ID)) 
  left_join(final_ref_rulengt) %>% 
  mutate(STR_length = STR_mean * RU.length) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_len사용
#write.table(c,"~/Desktop/KU/@research/STR/figure/figure1.notintersect_TRGT_EH_merge_processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")

c <- read.table("~/Desktop/KU/@research/STR/figure/figure1.notintersect_TRGT_EH_merge_processing.txt",header = T)
head(c)

c %>% ggplot(aes(x=g,y=STR_length,fill=g)) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  scale_fill_manual(values = c("TRGT" = "#F8766D","EH" = "#619CFF")) + 
  theme_step1()


#  이건 고민

c %>% ggplot(aes(x=g,y=GC)) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  theme_step1()


c %>% ggplot(aes(x=g,y=STR_length)) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  theme_step1()

head(merge_call_rate)
#c %>% 
#test
head(c)
c %>% mutate(GC_count = str_count(MOTIFS, "C") + str_count(MOTIFS, "G")) %>% 
  mutate(GC_count = GC_count*STR_mean) %>%
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>%
  ggplot(aes(x=g,y=GC_count)) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  labs(y="# of G or C") + 
  theme_step1()

  
  

c %>% filter(GC != 0) %>%
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>%
  ggplot(aes(x=g,y=GC*100)) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # 박스 플롯 오버레이
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") + # 평균값 추가
  labs(y="GC contents (%)") + 
  theme_step1()


c %>% filter(GC != 0) %>%
  ggplot(aes(x=g,y=GC)) + 
  geom_hex(bins=15)


c %>% ggplot(aes(x=g,y=GC)) + 
  geom_violin(trim = FALSE) +  # 바이올린 플롯
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # 박스 플롯 오버레이
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +  # 평균값 추가
  labs(x = "플랫폼", y = "GC 함량", title = "플랫폼별 GC 함량 밀도") +
  theme_minimal()

head(c)
c %>% #filter(STR_length < 100) %>% 
  ggplot(aes(x=GC,fill=g)) +
  #geom_histogram(position="identity", alpha=0.5)
  geom_bar(position = 'dodge',width = 0.05,)

c %>% filter(GC < 100) %>% 
  ggplot(aes(y=GC,fill=g)) +
  geom_histogram()

head(b)

head(tmp2)
#tmp2 %>% 
head(c)
  #head()
tmp2 %>% filter(is.na(g)) %>% head(10)
tmp2 %>% filter(STR_ID == "chr1_181015_181027")
tmp2 %>% count(g)
tmp2 %>% #filter(g == "EH") %>% head()
  select(STR_ID, MOTIFS,GC, EH_pass,EH_pass_range) %>% unique() -> c2
c2
venn_rawdata %>% filter(STR_DB == "chr10_100315009_100315119")
c2 %>% filter(is.na(EH_pass))
c %>% filter(g == "TRGT") %>% 
  left_join(c2) %>% head()
  
c %>% filter(g == "TRGT") %>% left_join(c2) %>% #head()
  ggplot(aes(x=EH_pass,y=GC,fill=EH_pass_range)) +
  geom_point()
#library(ggplot2)

c %>% filter(g == "TRGT") %>% left_join(c2) %>% #head()
  ggscatter(x="EH_pass",y= "GC",
            color = "black", shape = 21, size = 3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))

c %>% filter(g == "TRGT") %>% left_join(c2) %>% 
  ggplot(aes(x=EH_pass_range,y=GC)) +
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # 박스 플롯 오버레이
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red")  # 평균값 추가

head(c)
head(c2)

c2 %>% filter(STR_ID %in% c$STR_ID) #%>% na.omit()

c %>% filter(STR_ID %in% c2$STR_ID)
## 이거 굿굿
c %>% filter(g == "TRGT") %>% #head()
  left_join(c2) %>% #head()
  ggplot(aes(x=EH_pass,y=GC)) +
  geom_hex(bins = 15) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "STR call rate for SRS", y = "GC contents (%)", fill = "# of STRs") +  # 축 및 범례 라벨
  theme_step1()

#c %>% filter(g == "TRGT") %>% #select(STR_ID) %>%
 # left_join(c2) %>% write.table("~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect_withEHcallrate.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#rm(df)
trgt_withEHcall <- read_table("~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect_withEHcallrate.txt")
trgt_withEHcall

head(trgt_withEHcall)

trgt_withEHcall %>%
  ggplot(aes(x=EH_pass,y=GC)) +
  geom_hex(bins = 15) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "STR call rate for SRS", y = "GC contents (%)", fill = "# of STRs") +  # 축 및 범례 라벨
  theme_step1()

library(ggplot2)
#ggplot2::ggpl
c %>% filter(g == "TRGT") %>% #select(STR_ID) %>%
  left_join(c2) %>% head()
  ggplot(aes(x=EH_pass,y=GC)) +
  geom_density_2d_filled(alpha = 0.5) + 
  geom_density_2d(linewidth = 0.25, colour = "black") + 
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "STR call rate for SRS", y = "GC contents (%)") +  # 축 및 범례 라벨
  theme_step1()

c %>% mutate(GC_count = str_count(MOTIFS, "C") + str_count(MOTIFS, "G")) %>% 
  mutate(GC_count = GC_count*STR_mean) %>%
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>% #head()
  ggplot(aes(x=STR_length,y=GC)) + 
  geom_density_2d_filled(alpha = 0.5) + 
  geom_density_2d(linewidth = 0.25, colour = "black") + 
  labs(x = "STR length (bp)", y = "GC contents (%)") +  # 축 및 범례 라벨
  facet_grid(~g)
  

head(c)
c %>% mutate(GC_count = str_count(MOTIFS, "C") + str_count(MOTIFS, "G")) %>% 
  mutate(GC_count = GC_count*STR_mean) %>%
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>% #head()
  ggplot(aes(x=RU.length,y=GC)) + 
  geom_density_2d_filled(alpha = 0.5) + 
  geom_density_2d(linewidth = 0.25, colour = "black") + 
  labs(x = "RU.length (bp)", y = "GC contents (%)") +  # 축 및 범례 라벨
  facet_grid(~g)

head(c)
c %>% group_by(g) %>% summarise(max(STR_length))

c %>% mutate(GC_count = str_count(MOTIFS, "C") + str_count(MOTIFS, "G")) %>% 
  mutate(GC_count = GC_count*STR_mean) %>% #head()
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>% filter(STR_length <= 100) %>% 
  count(g)
  

c %>% mutate(GC_count = str_count(MOTIFS, "C") + str_count(MOTIFS, "G")) %>% 
  mutate(GC_count = GC_count*STR_mean) %>% #head()
  filter(GC != 0) %>%# head()
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>% #filter(GC_count >= 3) %>%
  #filter(STR_length <= 100) %>%
  mutate(STR_length = ifelse(STR_length >= 100,100,STR_length)) %>% #head()
  ggplot(aes(x=STR_length,y=GC*100)) + 
  geom_density2d() + 
  #geom_density_2d_filled(alpha = 0.5) + 
  geom_density_2d(linewidth = 0.25, colour = "black") + 
  labs(x = "STR length (bp)", y = "GC contents (%)") +  # 축 및 범례 라벨
  facet_grid(~g,scales = "free")

c %>% mutate(GC_count = str_count(MOTIFS, "C") + str_count(MOTIFS, "G")) %>% 
  mutate(GC_count = GC_count*STR_mean) %>% #head()
  filter(GC != 0) %>%# head()
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>% #filter(GC_count >= 3) %>%
  #filter(STR_length <= 100) %>%
  #mutate(STR_length = ifelse(STR_length >= 100,100,STR_length)) %>% #head()
  ggplot(aes(y=GC*100,x=g,fill=g)) + 
  geom_violin() + 
  #geom_density_2d_filled(alpha = 0.5) + 
  #geom_density_2d(linewidth = 0.25, colour = "black") + 
  labs(x = "STR length (bp)", y = "GC contents (%)") +  # 축 및 범례 라벨
  facet_grid(~g,scales = "free")


  
c %>% mutate(GC_count = str_count(MOTIFS, "C") + str_count(MOTIFS, "G")) %>% 
  mutate(GC_count = GC_count*STR_mean) %>% #head()
  filter(GC != 0) %>%# head()
  mutate(g = ifelse(g == "EH","SRS","LRS")) %>% #filter(GC_count >= 50) %>%
  #filter(STR_length <= 100) %>%
  mutate(STR_length = ifelse(STR_length >= 50,50,STR_length)) %>% 
  ggplot(aes(x=STR_length,fill=g)) + 
  geom_bar(position = 'dodge',width = 0.05)





head(c2)
c %>% filter(g == "TRGT") %>% left_join(c2) %>% #head()
  ggplot(aes(x=EH_pass_range,y=GC)) +
  geom_hex(bins = 15) +
  scale_fill_gradient(low = "blue", high = "red")
###3
head(c)
c %>% left_join(merge_call_rate)
head(final_ref_pro)


### pandepth GC rerererererer
notintersect_TRGT_EH_merge_processing<- read.table("~/Desktop/KU/@research/STR/figure/figure1.notintersect_TRGT_EH_merge_processing.txt",header = T)
notintersect_TRGT_EH_merge_processing %>% filter(g=="EH") %>% select(STR_ID) -> eh_not
notintersect_TRGT_EH_merge_processing %>% filter(g!="EH") %>% select(STR_ID) -> trgt_not

head(notintersect_TRGT_EH_merge_processing)
head(c)
eh_notinter_pandepth <- read_table("~/Desktop/KU/@research/STR/eh/pandepth.merge.from.eh.oribam.txt") %>% filter(STR_ID %in% eh_not$STR_ID)
trgt_notinter_pandepth <- read_table("~/Desktop/KU/@research/STR/trgt/pandepth/pandepth.merge.from.trgt.oribam.txt") %>% filter(STR_ID %in% trgt_not$STR_ID)
head(eh_notinter_pandepth)
head(trgt_notinter_pandepth)
trgt_notinter_pandepth %>% select(STR_ID,`GC(%)`,MeanDepth,ID) -> trgt_notinter_pandepth_a
eh_notinter_pandepth %>% select(STR_ID,`GC(%)`,MeanDepth,ID) -> eh_notinter_pandepth_a

#eh %>% filter(STR_ID %in% not_intersect$STR_DB) %>% filter(ID != "NIH20N2594890") -> eh_notintersect
#trgt %>% filter(TRID %in% not_intersect$STR_DB) %>% filter(ID != "NIH23J3904558") -> trgt_notintersect
head(trgt_notinter)
head(trgt_notinter_pandepth_a)
trgt_notinter <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect.txt") %>% rename(STR_ID = TRID) 
trgt_notinter %>% left_join(trgt_notinter_pandepth_a) %>% filter(ID != "NIH23J3904558") %>% filter(!is.na(`GC(%)`)) %>%
  select(ID,STR_ID,`GC(%)`,MeanDepth) -> trgt_notinter_pandepth_b


eh_notinter <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure1.eh_notintersect_processing.txt") 
#head(eh_notinter)

trgt_notinter %>% left_join(eh_notinter_pandepth_a) %>% filter(ID != "NIH23J3904558") %>% filter(!is.na(`GC(%)`)) %>%
  select(ID,STR_ID,`GC(%)`,MeanDepth) -> trgt_notinter_pandepth_b


library(ComplexHeatmap)
  
### annotation annova

anno <- read_table("~/Desktop/KU/@research/STR/db/annovar/Raw.anno.merge.processing.onlyneed.txt")
trgt_withEHcall <- read_table("~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect_withEHcallrate.txt")
c <- read.table("~/Desktop/KU/@research/STR/figure/figure1.notintersect_TRGT_EH_merge_processing.txt",header = T)

#chr2_140766822_140766844
(140766822 + 140766844)/2
(100001401 + 100001413)/2
head(trgt_withEHcall)
head(anno)
trgt_withEHcall %>% left_join(anno) %>% count(type) #%>%

trgt_withEHcall %>% left_join(anno) %>% 
  ggplot(aes(x=type,y=EH_pass)) + 
  geom_violin()

trgt_withEHcall %>% left_join(anno) %>% 
  ggplot(aes(x=main_cyto_type,y=EH_pass)) + 
  geom_violin()


  #count(main_cyto_type)
############ annotation vep

annot_patho <- read_table("~/Desktop/KU/@research/STR/db/annotation/catalog.GRCh38.with_adjacent_repeats.TRGT.vep",skip = 40)
annot <- read_table("~/Desktop/KU/@research/STR/db/annotation/eh.v5_w_gangstr.v13.target.region.sort.vep",skip = 40)

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)


head(annot_patho)
head(annot)

annot %>% select(ID,Consequence)

colnames(annot_patho)[1] <- "ID"
colnames(annot)[1] <- "ID"

head(annot)
annot %>% select(ID,Consequence) 
annot %>% select(ID,Consequence) %>% filter(str_detect(ID,"chr")) %>%
  filter(ID %in%final_ref$ID) -> annot_v2
head(annot_v2)
annot_v2 %>% count(Feature_type)

annot_patho %>% select(ID,Consequence) -> annot_patho2

annot_v2
annot_v2 %>% count(ID,Consequence) %>% filter(n != 1)
annot_patho %>% select(ID,Location,Feature_type,Consequence) %>% unique() %>%
  rbind(annot %>% select(ID,Location,Feature_type,Consequence) %>% unique()) -> anno_merge

annot_patho2 %>%
  tidyr::separate_rows(Consequence, sep = ",") -> annot_patho3
head(annot_patho2)

annot_v2 %>%
  tidyr::separate_rows(Consequence, sep = ",") -> annot_v3

head(annot_v3)
annot_patho2  %>% select(Consequence) %>% unique() -> annot_patho4
annot_v3 %>% select(Consequence) %>% unique() -> anno_seq
head(anno_seq)
head(annot_patho2)

annot_patho3 %>% rbind(anno_seq) %>% unique()


data_frame(
  Consequence  = c("coding_sequence_variant","3_prime_UTR_variant","5_prime_UTR_variant","downstream_gene_variant",
  "upstream_gene_variant","intron_variant","intergenic_variant","start_lost","non_coding_transcript_exon_variant","non_coding_transcript_variant",
  "NMD_transcript_variant",
  "mature_miRNA_variant","-"),
type = c("Exon","3-UTR","5-UTR","downstream","upstream","intron","intergenic","start_log","other","other","other","other",'other'),
rank = c(seq(1:13))
) -> priority
head(annot_patho2)
head(annot_patho2)
head(annot_patho2)
annot_patho3 %>% unique %>% left_join(priority) %>% group_by(ID) %>%
  filter(rank == min(rank)) -> annot_patho_qc
  

annot_v3 %>% unique %>% left_join(priority) %>% group_by(ID) %>%
  filter(rank == min(rank)) -> annot_v3_qc

head(annot_v3_qc)
head(annot_patho_qc)
head(onlyEH)
library(ggpie)
library(ggplot2)
head(venn_rawdata)

venn_rawdata %>% na.omit() -> common_STR
final_ref %>% filter(ID %in% venn_rawdata$TRGT) %>% filter(!(ID %in% venn_rawdata$EH)) %>% select(ID) -> onlyTRGT_ID
final_ref %>% filter(ID %in% venn_rawdata$EH) %>% filter(!(ID %in% venn_rawdata$TRGT)) %>% select(ID) -> onlyEH_ID


annot_v3_qc %>% rbind(annot_patho_qc) %>% filter(ID %in% common_STR$STR_DB) -> annot_common_str
annot_v3_qc %>% rbind(annot_patho_qc) %>% filter(ID %in% onlyTRGT_ID$ID) -> annot_onlyTRGT_str
annot_v3_qc %>% rbind(annot_patho_qc) %>% filter(ID %in% onlyEH_ID$ID) -> annot_onlyEH_str

head(annot_common_str)
annot_common_str %>% select(-rank) %>% write.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.CommonSTR.forpie.txt",row.names = F,col.names = T,quote = F,sep = "\t")
annot_onlyTRGT_str %>% select(-rank) %>% write.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.onlyTRGT_STR.forpie.txt",row.names = F,col.names = T,quote = F,sep = "\t")
annot_onlyEH_str %>% select(-rank) %>% write.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.onlyEH_STR.forpie.txt",row.names = F,col.names = T,quote = F,sep = "\t")


ggdonut(annot_common_str,group_key = 'type',count_type = 'full', 
        label_info = c("group","ratio"), label_type = "horizon",
        label_size = 4, label_pos = "in", label_threshold = 10) + theme(legend.position = 'none')

######
