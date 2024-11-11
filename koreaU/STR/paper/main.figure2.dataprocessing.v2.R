### main figure 2 data processing
library(tidyverse)
library(data.table)
library(ggbreak)
library(ggpubr)
library(ggplot2)
library(cowplot)

library(grid)
library(ggpie)
library(ggforce)

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro


head(final_ref)
final_ref %>% filter(!str_detect(ID,"chr"))

#save(eh_simple_pass_intersect, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_simple_pass_intersect.RData")
#save(eh_simple_pass_intersect, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_simple_pass_intersect.RData")
#save(eh_complex_pass, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_complex_pass.RData")

#save(trgt_simple_pass_intersect, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_simple_pass_intersect.RData")
#save(trgt_simple_pass_onlytrgt, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_simple_pass_onlytrgt.RData")
#save(trgt_complex_pass, file = "/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_complex_pass.RData")
load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_complex_pass.RData")
head(trgt_complex_pass[1:5])
trgt_complex_pass %>% select(TRID) %>% unique() %>% dim()
trgt_complex_pass %>% select(TRID) %>% unique()

trgt_complex_pass %>% filter(ID == "NIH23J3904558")
eh_simple_pass_intersect_gt %>% filter(ID == "NIH20N2594890")
head(final_ref)
final_ref %>% filter(ID %in% trgt_complex_pass$TRID)
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR


load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_simple_pass_intersect.RData")
load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_simple_pass_intersect.RData")

head(eh_simple_pass_intersect)
head(trgt_simple_pass_intersect)

table(eh_simple_pass_intersect$FILTER)

eh_simple_pass_intersect %>% count(STR_ID %in% common_STR$STR_DB)

common_STR %>% count(STR_DB %in% eh_simple_pass_intersect$STR_ID)



head(eh_simple_pass_intersect_gt)
eh_simple_pass_intersect %>% select(-type,-ADSP,-ADFL,-ADIR) %>% #filter(FILTER == "PASS") %>%
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>%
  mutate(ALT = ifelse(ALT == ".",0,ALT)) %>% mutate(GT = str_split_fixed(GT,"/",2),ALT= str_split_fixed(ALT,",",2)) -> eh_simple_pass_intersect_gt

eh_simple_pass_intersect_gt %>% mutate(STR1 = ifelse(GT[,1] == "1",ALT[,1],ifelse(GT[,1] == "2",ALT[,2],REF))) %>%
  mutate(STR2 = ifelse(GT[,2] == "1",ALT[,1],ifelse(GT[,2] == "2",ALT[,2],REF))) %>% select(-ALT,-GT,-REF) %>% mutate(STR1 = as.integer(STR1),STR2 = as.integer(STR2)) -> eh_simple_pass_intersect_gt

head(eh_simple_pass_intersect)

eh_simple_pass_intersect %>% select(ID,STR_ID,type) %>% mutate(type1 = str_split_fixed(type,"/",2)[,1],type2 = str_split_fixed(type,"/",2)[,2]) %>% 
  select(-type) -> eh_simple_pass_intersect_type

head(eh_simple_pass_intersect_gt)
head(eh_simple_pass_intersect_type)

eh_simple_pass_intersect_gt %>% select(-FILTER) %>% left_join(eh_simple_pass_intersect_type) -> eh_simple_pass_intersect_processing

eh_simple_pass_intersect_processing %>% mutate(EH_STR1 = ifelse(STR1 < STR2, STR1,STR2),EH_STR2 = ifelse(STR1 < STR2, STR2,STR1)) %>%
  mutate(EH_type1 = ifelse(STR1 < STR2, type1,type2),EH_type2 = ifelse(STR1 < STR2, type2,type1)) %>% 
  select(-STR1,-STR2,-type1,-type2)-> eh_simple_pass_intersect_processing

head(eh_simple_pass_intersect_processing)
eh_simple_pass_intersect_processing %>% filter(ID != "NIH20N2594890") -> eh_simple_pass_intersect_processing

eh_simple_pass_intersect_processing %>% select(-CHROM,-POS) %>% write.table("~/Desktop/KU/@research/STR/figure/figure2/eh_simple_pass_intersect_processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")



###
head(trgt_simple_pass_intersect)

trgt_simple_pass_intersect %>% select(ID,TRID,MOTIFS,Allele,MC) %>% #head()
  pivot_wider(names_from = Allele,values_from = MC) %>% rename("MC1" = allele_1,"MC2" = allele_2) -> trgt_simple_pass_intersect_gt

trgt_simple_pass_intersect %>% select(ID,TRID,MOTIFS,Allele,AP) %>% #head()
  pivot_wider(names_from = Allele,values_from = AP) %>% rename("AP1" = allele_1,"AP2" = allele_2) -> trgt_simple_pass_intersect_AP

trgt_simple_pass_intersect %>% select(ID,TRID,MOTIFS,Allele,AM) %>% #head()
  pivot_wider(names_from = Allele,values_from = AM) %>% rename("AM1" = allele_1,"AM2" = allele_2) -> trgt_simple_pass_intersect_AM


head(trgt_simple_pass_intersect_gt)
head(trgt_simple_pass_intersect_AM)

trgt_simple_pass_intersect_gt %>% left_join(trgt_simple_pass_intersect_AM) %>% left_join(trgt_simple_pass_intersect_AP) -> trgt_simple_pass_intersect_merge

head(trgt_simple_pass_intersect_merge)
head(trgt_simple_pass_intersect_merge)
trgt_simple_pass_intersect_merge %>% filter(MC1 == ".")
trgt_simple_pass_intersect_merge %>% filter(MC2 == ".")
trgt_simple_pass_intersect_merge %>% filter(is.na(AP1))
trgt_simple_pass_intersect_merge %>% filter(is.na(AP2))
trgt_simple_pass_intersect_merge %>% filter(AP2 == ".")
trgt_simple_pass_intersect_merge %>% filter(AP1 == ".")
trgt_simple_pass_intersect_merge %>% filter(AP1 == "0")

trgt_simple_pass_intersect_merge %>%
  mutate(across(AM1:AP2, ~ ifelse(. == ".", -1, .))) ->  trgt_simple_pass_intersect_merge
head(trgt_simple_pass_intersect_merge)

trgt_simple_pass_intersect_merge %>% mutate(MC1 = as.numeric(MC1),MC2 = as.numeric(MC2)) %>%
  mutate(TRGT_STR1 = ifelse(MC1 < MC2, MC1,MC2),TRGT_STR2 = ifelse(MC1 < MC2, MC2,MC1)) %>%
  mutate(TRGT_AM1 = ifelse(MC1 < MC2, AM1,AM2),TRGT_AM2 = ifelse(MC1 < MC2, AM2,AM1)) %>%
  mutate(TRGT_AP1 = ifelse(MC1 < MC2, AP1,AP2),TRGT_AP2 = ifelse(MC1 < MC2, AP2,AP1)) %>% 
  select(-MC1,-MC2,-AM1,-AM2,-AP1,-AP2)  -> trgt_simple_pass_intersect_merge_processing

head(trgt_simple_pass_intersect_merge_processing)
trgt_simple_pass_intersect_merge_processing %>% count(TRGT_STR1 < TRGT_STR2)
trgt_simple_pass_intersect_merge_processing %>% count(TRGT_STR1 <= TRGT_STR2)

trgt_simple_pass_intersect_merge_processing %>% filter(ID != "NIH23J3904558") -> trgt_simple_pass_intersect_merge_processing
head(trgt_simple_pass_intersect_merge_processing)

#write.table(trgt_simple_pass_intersect_merge_processing,"~/Desktop/KU/@research/STR/figure/figure2/trgt_simple_pass_intersect_merge_processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")

######################################
trgt_simple_pass_intersect_merge_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure2/trgt_simple_pass_intersect_merge_processing.txt")
eh_simple_pass_intersect_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_simple_pass_intersect_processing.txt")
head(eh_simple_pass_intersect_processing)
trgt_simple_pass_intersect_merge_processing %>% count(TRGT_STR1<=TRGT_STR2)
eh_simple_pass_intersect_processing %>% count(EH_STR1<=EH_STR2)
trgt_simple_pass_intersect_merge_processing %>% filter(!(TRGT_STR1<=TRGT_STR2))

head(trgt_simple_pass_intersect_merge_processing)
head(eh_simple_pass_intersect_processing)
trgt_simple_pass_intersect_merge_processing %>% select(ID) %>% unique() %>% dim()
eh_simple_pass_intersect_processing %>% select(ID) %>% unique() %>% dim()

trgt_simple_pass_intersect_merge_processing %>% select(TRID) %>% unique() %>% dim()
eh_simple_pass_intersect_processing %>% select(ID) %>% unique() %>% dim()


ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)

ref %>% select(Revio,Illumina) -> ref
head(ref)
colnames(ref) <- c("ID","EH_ID")

eh_simple_pass_intersect_processing  %>% rename(EH_ID = ID) %>% left_join(ref) %>%
  select(STR_ID,RU,ID,EH_STR1,EH_STR2) -> eh_simple_pass_intersect_processing_qt

head(eh_simple_pass_intersect_processing_qt)
trgt_simple_pass_intersect_merge_processing %>% select(-MOTIFS) %>% rename(STR_ID = TRID) %>% left_join(eh_simple_pass_intersect_processing_qt) -> a
head(a)  

#write.table(a,"~/Desktop/KU/@research/STR/figure/figure2/eh_trgt_merge_simple_pass_intersect_forConcordance.txt",col.names = T,row.names = F,quote = F,sep = "\t")

head(trgt_simple_pass_intersect_merge_processing)

trgt_simple_pass_intersect_merge_processing %>% select(ID,TRID,TRGT_AP1,TRGT_AP2) %>% #filter(TRGT_AP1 == 0) %>% 
  pivot_longer(3:4) %>% #filter(value >= 0.5) %>% head()
  mutate(allele = ifelse(name == "TRGT_AP1",1,2)) %>% select(-name)-> a
head(a)
colnames(a) <- c("ID","STR_ID","AP","allele") 
write.table(a[,c("ID","STR_ID","allele","AP")],"~/Desktop/KU/@research/STR/trgt/trgt_simpleSTR_intersect.AP.score.byAllele.txt",col.names = T,row.names = F,quote = F,sep = "\t")

a %>% rename(STR_ID = TRID) %>% left_join(concordance_rate) %>% 
  mutate(g = ifelse(APmean >=0.5,"0.5>=","<0.5")) %>% 
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "0",
    concordance_rate > 0 & concordance_rate < 0.1 ~ "(0,0.1)",
    concordance_rate >= 0.1 & concordance_rate < 0.2 ~ "[0.1,0.2)",
    concordance_rate >= 0.2 & concordance_rate < 0.3 ~ "[0.2,0.3)",
    concordance_rate >= 0.3 & concordance_rate < 0.4 ~ "[0.3,0.4)",
    concordance_rate >= 0.4 & concordance_rate < 0.5 ~ "[0.4,0.5)",
    concordance_rate >= 0.5 & concordance_rate < 0.6 ~ "[0.5,0.6)",
    concordance_rate >= 0.6 & concordance_rate < 0.7 ~ "[0.6,0.7)",
    concordance_rate >= 0.7 & concordance_rate < 0.8 ~ "[0.7,0.8)",
    concordance_rate >= 0.8 & concordance_rate < 0.9 ~ "[0.8,0.9)",
    concordance_rate >= 0.9 & concordance_rate < 1 ~ "[0.9,1)",
    concordance_rate == 1 ~ "1")) %>% count(concordance_rate_range,g) %>% 
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                       "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                       "[0.8,0.9)", "[0.9,1)","1")),y=n,fill=g)) + 
  geom_bar(stat = "identity",position = 'fill')
##### concordace

eh_trgt_merge_simple_pass_intersect_forConcordance <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_trgt_merge_simple_pass_intersect_forConcordance.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance)
eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(STR_ID) %>% unique() %>% dim()

eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID:TRGT_STR2,RU:EH_STR2) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(concordance_point = STR1 + STR2) %>% group_by(STR_ID) %>% 
  summarise(concordance_point_sum = sum(concordance_point)) %>% 
  mutate(concordance_rate = concordance_point_sum/130) -> concordance_bySTR_simpleSTR
  
concordance_bySTR_simpleSTR
concordance_bySTR_simpleSTR %>% count(str_detect(STR_ID,"chr"))

# general 311470 
# patho 63


concordance_bySTR_simpleSTR %>% filter(concordance_rate == 1) %>%
  count(str_detect(STR_ID,"chr"))

# general 170185 
# patho 18



####
head(eh_trgt_merge_simple_pass_intersect_forConcordance)





concordance_bySTR_simpleSTR

concordance_bySTR_simpleSTR %>% mutate(STR_type = ifelse(str_detect(STR_ID,"chr"),"Normal","Pathogenic")) %>%
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "0",
    concordance_rate > 0 & concordance_rate < 0.1 ~ "(0,0.1)",
    concordance_rate >= 0.1 & concordance_rate < 0.2 ~ "[0.1,0.2)",
    concordance_rate >= 0.2 & concordance_rate < 0.3 ~ "[0.2,0.3)",
    concordance_rate >= 0.3 & concordance_rate < 0.4 ~ "[0.3,0.4)",
    concordance_rate >= 0.4 & concordance_rate < 0.5 ~ "[0.4,0.5)",
    concordance_rate >= 0.5 & concordance_rate < 0.6 ~ "[0.5,0.6)",
    concordance_rate >= 0.6 & concordance_rate < 0.7 ~ "[0.6,0.7)",
    concordance_rate >= 0.7 & concordance_rate < 0.8 ~ "[0.7,0.8)",
    concordance_rate >= 0.8 & concordance_rate < 0.9 ~ "[0.8,0.9)",
    concordance_rate >= 0.9 & concordance_rate < 1 ~ "[0.9,1)",
    concordance_rate == 1 ~ "1")) %>% #mutate(concordance_match = ifelse(concordance_rate == 1))
  #group_by(concordance_rate_range) %>%
  count(concordance_rate_range) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance.range.bySTR_simpleSTR.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  mutate(facet = ifelse(concordance_rate_range %in% c("[0.9,1)","1"),2,1)) %>% #head()
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                              "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                              "[0.8,0.9)", "[0.9,1)","1")),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Range of") + 
  geom_text(aes(label = scales::comma(n), y = n), # 바의 위에 텍스트를 배치하기 위해 y를 n보다 약간 크게 설정
            size = 6, family = "Arial", vjust = 0) +
  facet_row(~ facet, scales = "free", space = "free") + 
    theme(strip.text = element_blank())
  
  
eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID:TRGT_STR2,RU:EH_STR2) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(concordance_point = (STR1 + STR2)) %>% group_by(ID) %>% 
  summarise(concordance_point_sum = sum(concordance_point)) %>% 
  mutate(concordance_rate = concordance_point_sum/(311542*2)) -> concordance_byID 

  
head(concordance_byID)

#concordance_byID %>% write.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance.range.byID_simpleSTR.txt",col.names = T,row.names = F,quote = F,sep = "\t")

concordance_byID %>% ggplot(aes(x = "", y = concordance_rate)) + 
  geom_violin(trim=FALSE,fill="gray") + 
  labs(y="# of STRs") +
  theme_step1() + 
  theme(axis.ticks = element_blank(),
        axis.title.x = element_blank())

ru.scale = c(seq(2,19),"20+")
ru.scale1 = c(seq(2,6))
ru.scale2 = c(seq(7,19),"20+")


head(concordance_bySTR_simpleSTR) 
head(final_ref)
concordance_bySTR_simpleSTR %>% left_join(final_ref_pro) %>%  #head()
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "x = 0",
    concordance_rate > 0 & concordance_rate < 1 ~ "0 < x < 1",
    concordance_rate == 1 ~ "x = 1")) %>%
  count(RU.length,concordance_rate_range) %>%
  #mutate(facet = ifelse(RU.length %in% ru.scale1,1,2)) %>% #head() 
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n,fill=factor(concordance_rate_range,levels=c("x = 0","0 < x < 1","x = 1"))))  +
  geom_bar(stat='identity',position = 'fill') + 
  labs(x = "Length of Repeat Units", y = "Proportion", fill = "Concordance rate") +
  theme_step1() + 
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
  

concordance_bySTR_simpleSTR %>% left_join(final_ref_pro) %>%  #head()
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "x = 0",
    concordance_rate > 0 & concordance_rate < 1 ~ "0 < x < 1",
    concordance_rate == 1 ~ "x = 1")) %>%
  count(RU.length,concordance_rate_range) %>% #head()
  group_by(RU.length) %>%
  mutate(prop = prop.table(n)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=prop,fill=factor(concordance_rate_range,levels=c("x = 0","0 < x < 1","x = 1"))))  +
  geom_bar(stat='identity') +
  labs(x = "Length of Repeat Units", y = "Proportion", fill = "Concordance rate") +
  theme_step1() + 
  coord_cartesian(ylim = c(0.98, 1)) +   # y축을 0.99~1로 제한 
  theme(legend.position = "none",
        axis.title.x = element_blank())

concordance_bySTR_simpleSTR %>% left_join(final_ref_pro) %>%  #head()
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "x = 0",
    concordance_rate > 0 & concordance_rate < 1 ~ "0 < x < 1",
    concordance_rate == 1 ~ "x = 1")) %>%
  count(RU.length,concordance_rate_range) %>% #head()
  group_by(RU.length) %>%
  mutate(prop = prop.table(n)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=prop,fill=factor(concordance_rate_range,levels=c("x = 0","0 < x < 1","x = 1"))))  +
  geom_bar(stat='identity') +
  labs(x = "Length of Repeat Units", y = "Proportion", fill = "Concordance rate") +
  theme_step1() + 
  coord_cartesian(ylim = c(0.98, 1)) +   # y축을 0.99~1로 제한 
  theme(legend.position = "bottom",
        legend.direction = "horizontal")


  
  
ru.scale = c(seq(2,19),"20+")
ru.scale1 = c(seq(2,6))
ru.scale2 = c(seq(7,19),"20+")

library(ggplot2)
library(patchwork)
head(concordance_bySTR_simpleSTR)
#write.table(concordance_bySTR_simpleSTR,"~/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt",col.names = T,row.names = F,quote = F,sep = "\t")
concordance_bySTR_simpleSTR <- read.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt",header = T)
# 0~0.99 구간을 그리는 그래프
concordance_bySTR_simpleSTR %>% 
  left_join(final_ref_pro) %>%
  mutate(RU.length = ifelse(RU.length %in% ru.scale, RU.length, "20+")) %>%
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "x = 0",
    concordance_rate > 0 & concordance_rate < 1 ~ "0 < x < 1",
    concordance_rate == 1 ~ "x = 1"
  )) %>%
  count(RU.length, concordance_rate_range) %>%
  group_by(RU.length) %>%
  mutate(prop = prop.table(n)) %>%
  ggplot(aes(x = factor(RU.length, levels = ru.scale), y = prop, fill = factor(concordance_rate_range, levels = c("x = 0", "0 < x < 1", "x = 1")))) +
  geom_bar(stat = 'identity') +
  labs(x = "Length of Repeat Units", y = "Proportion", fill = "Concordance rate") +
  geom_text(aes(label = ifelse(concordance_rate_range %in% c("x = 0"),n,""),y = 1,hjust = 0),size = 4) + 
  theme_step1() +
  theme(legend.direction = "horizontal",
    legend.position = "bottom") +
  coord_flip()-> p1

p1
# 0.99~1 구간을 확대해서 그리는 그래프
p2 <- concordance_bySTR_simpleSTR %>% 
  left_join(final_ref_pro) %>%
  mutate(RU.length = ifelse(RU.length %in% ru.scale, RU.length, "20+")) %>%
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "x = 0",
    concordance_rate > 0 & concordance_rate < 1 ~ "0 < x < 1",
    concordance_rate == 1 ~ "x = 1"
  )) %>%
  count(RU.length, concordance_rate_range) %>%
  group_by(RU.length) %>%
  mutate(prop = prop.table(n)) %>%
  ggplot(aes(x = factor(RU.length, levels = ru.scale), y = prop, fill = factor(concordance_rate_range, levels = c("x = 0", "0 < x < 1", "x = 1")))) +
  geom_bar(stat = 'identity') +
  labs(x = "~", y = "Proportion", fill = "Concordance rate") +
  theme_step1() +
  coord_cartesian(ylim = c(0.98, 1)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.01))  + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y= element_blank(),
        axis.line.x = element_blank())

# 두 그래프를 위아래로 결합
cowplot::plot_grid(p2,p1,nrow = 2,rel_heights = c(1,1.5))


(p2 / p1) + 
  ylab("111")  # y축에 글자를 추가######### concordance 특징

head(final_ref)
final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro

eh_trgt_merge_simple_pass_intersect_forConcordance <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_trgt_merge_simple_pass_intersect_forConcordance.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance)


eh_trgt_merge_simple_pass_intersect_forConcordance %>% #select(ID:TRGT_STR2,RU:EH_STR2) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% rename(MOTIFS = RU) %>%
  filter(STR1 == 0, STR2 == 0) %>% select(ID:MOTIFS) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  #mutate(Whole.length = RU.length*mean_repeat_count) %>% #head()
  #group_by(STR,MOTIFS) %>%
  #summarise(RU.length = mean(RU.length),mean_Whole.length = mean(Whole.length),mean_repeat_count= mean(mean_repeat_count)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> eh_trgt_merge_simple_pass_intersect_forConcordance_match_info

head(eh_trgt_merge_simple_pass_intersect_forConcordance_match_info)


eh_trgt_merge_simple_pass_intersect_forConcordance_match_info %>% filter(TRGT_AP1 == -1)
  #mutate(mean_repeat_count = (TRGT_STR1+TRGT_STR1)/2) %>%

eh_trgt_merge_simple_pass_intersect_forConcordance_match_info %>%
  mutate(mean_repeat_count = (TRGT_STR1+TRGT_STR1)/2) %>% 
  mutate(Whole.length = RU.length*mean_repeat_count) %>% 
  group_by(STR_ID,MOTIFS) %>% 
  summarise(RU.length = mean(RU.length),mean_Whole.length = mean(Whole.length),mean_repeat_count= mean(mean_repeat_count)) %>% 
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> a


head(a)

a %>% ggplot(aes(x=mean_Whole.length,y=GC)) + 
  geom_point()

head(eh_trgt_merge_simple_pass_intersect_forConcordance)
eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID,STR_ID,TRGT_STR1,TRGT_STR2,EH_STR1,EH_STR2) %>% 
  pivot_longer(TRGT_STR1:EH_STR2) %>% mutate(name = ifelse(str_detect(name,"TRGT"),"TRGT","EH")) %>%
  group_by(STR_ID,name) %>%  #head()
  left_join(final_ref_pro) -> eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount

head(eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount)

eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount %>%
  mutate(Whole.length = RU.length*value) %>%
  summarise(mean_Whole.length = mean(Whole.length),# 중앙값
            q1_Whole.length = quantile(Whole.length, 0.25, na.rm = TRUE),# 1사분위수(1/4 값)
            q3_Whole.length = quantile(Whole.length, 0.75, na.rm = TRUE),
            max_Whole.length = max(Whole.length, na.rm = TRUE),         # 최대값
            min_Whole.length = min(Whole.length, na.rm = TRUE)) -> eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength


#write.table(eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength,"~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.txt",col.names = T,row.names = F,quote = F,sep = "\t")
eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength)


head(eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength)

eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% #ungroup() %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% head(1000) %>%
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = min_Whole.length_EH, ymax = max_Whole.length_EH), width = 0,alpha = 0.5) +  # y 방향 에러바
  geom_errorbarh(aes(xmin = min_Whole.length_TRGT, xmax = max_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "mean of STR length by LRS", y = "mean of STR length by SRS",title = "mean +max,-min") +
  theme_minimal() + 
  coord_fixed(ratio = 1) 


eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% #ungroup() %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% head(1000) %>%
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = min_Whole.length_EH, ymax = max_Whole.length_EH), width = 0,alpha = 0.5) +  # y 방향 에러바
  geom_errorbarh(aes(xmin = min_Whole.length_TRGT, xmax = max_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "median of STR length by LRS", y = "median of STR length by SRS",title = "median +max,-min") +
  theme_minimal() + 
  coord_fixed(ratio = 1) 

head(eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength)




eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% #ungroup() %>%
  pivot_wider(names_from = name,values_from = median_Whole.length:q3_Whole.length) %>% head(1000) %>%
  ggplot(aes(x = median_Whole.length_TRGT, y = median_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = q1_Whole.length_EH, ymax = q3_Whole.length_EH), width = 0,alpha = 0.5) +  # y 방향 에러바
  geom_errorbarh(aes(xmin = q1_Whole.length_TRGT, xmax = q3_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "Mean of STR length by LRS", y = "Mean of STR length by SRS") +
  theme_minimal() + 
  coord_fixed(ratio = 1) 



eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID,STR_ID,TRGT_STR1,TRGT_STR2,EH_STR1,EH_STR2) %>%
  pivot_longer(TRGT_STR1:EH_STR2) %>% mutate(name = ifelse(str_detect(name,"TRGT"),"TRGT","EH")) %>%
  left_join(final_ref_pro) -> eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs

head(eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs)
eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs %>% 
  mutate(Whole.length = RU.length * value) %>%
  group_by(name,MOTIFS) %>% 
  summarise(mean_Whole.length = mean(Whole.length, na.rm = TRUE),   # 중앙값
            q1_Whole.length = quantile(Whole.length, 0.25, na.rm = TRUE),# 1사분위수(1/4 값)
            q3_Whole.length = quantile(Whole.length, 0.75, na.rm = TRUE),
            max_Whole.length = max(Whole.length, na.rm = TRUE),         # 최대값
            min_Whole.length = min(Whole.length, na.rm = TRUE)) -> eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0
  
head(eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0)
#write.table(eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0,"~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.byMOTIFs.mean.txt",col.names = T,row.names = F,quote = F,sep = "\t")
eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.byMOTIFs.mean.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0)
eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% head(1000) %>% #head()
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = min_Whole.length_EH, ymax = max_Whole.length_EH), width = 0,러바
  geom_errorbarh(aes(xmin = min_Whole.length_TRGT, xmax = max_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "mean of STR length by LRS", y = "mean of STR length by SRS",title = "mean +max,-min") +
  theme_minimal() + 
  coord_fixed(ratio = 1) 


eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% #head(1000) %>% #head()
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  #geom_errorbar(aes(ymin = min_Whole.length_EH, ymax = max_Whole.length_EH), width = 0,alpha = 0.5) +  # y 방향 에러바
  #geom_errorbarh(aes(xmin = min_Whole.length_TRGT, xmax = max_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "mean of STR length by LRS", y = "mean of STR length by SRS",title = "mean by MOTIFs") +
  theme_minimal() + 
  coord_fixed(ratio = 1) 


eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% #head(1000) %>% #head()
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = min_Whole.length_EH, ymax = max_Whole.length_EH), width = 0,alpha = 0.5) +  # y 방향 에러바
  geom_errorbarh(aes(xmin = min_Whole.length_TRGT, xmax = max_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "mean of STR length by LRS", y = "mean of STR length by SRS",title = "mean by MOTIFs + max/min") +
  theme_minimal() + 
  coord_fixed(ratio = 1) 

head(eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0)
eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% #head(1000) %>% #head()
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = q1_Whole.length_EH, ymax = q3_Whole.length_EH), width = 0,alpha = 0.5) +  # y 방향 에러바
  geom_errorbarh(aes(xmin = q1_Whole.length_TRGT, xmax = q3_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "mean of STR length by LRS", y = "mean of STR length by SRS",title = "mean by MOTIFs + q1/q3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  theme_minimal() + 
  coord_fixed(ratio = 1)


#######
eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.txt")
eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.byMOTIFs.mean.txt")

head(eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength)
head(eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0)
simple_concordance <- read.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt",header = T)

simple_concordance %>% mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>% group_by(g) %>%
  summarise(mean(concordance_rate))
head(simple_concordance)
head(final_ref_pro)
simple_concordance %>% left_join(final_ref_pro) %>% group_by(MOTIFS) %>% summarise(mean_concordacne = mean(concordance_rate)) %>%
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> a

head(a)

eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 %>% left_join(a) %>% #head()
  filter(name == "TRGT") %>% 
  filter(mean_concordacne != 1) %>%  #head(1000) %>% #head()
  ggplot(aes(x=mean_concordacne,y=mean_Whole.length)) + 
  geom_hex(bins = 15) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "Concordance rate", y = "mean of STR length", fill = "# of STRs") +  # 축 및 범례 라벨
  theme_step1()

  
eh_trgt_merge_simple_pass_intersect_forConcordance_meanRepeatcount_byMOTIFs0 %>% left_join(a) %>% #head()
  filter(name == "TRGT") %>% 
  #filter(mean_concordacne != 1) %>%  #head(1000) %>% #head()
  ggplot(aes(x=mean_concordacne,y=GC)) + 
  geom_hex(bins = 15) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "Concordance rate", y = "mean of STR length", fill = "# of STRs") +  # 축 및 범례 라벨
  theme_step1()

eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% left_join(simple_concordance) %>% #head()
  filter(concordance_rate != 1) %>% select(-concordance_point_sum,-concordance_rate) %>% #head(1000) %>%head()
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% #head(1000) %>% #head()
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = q1_Whole.length_EH, ymax = q3_Whole.length_EH), width = 0,alpha = 0.5) +  # y 방향 에러바
  geom_errorbarh(aes(xmin = q1_Whole.length_TRGT, xmax = q3_Whole.length_TRGT), height = 0,alpha = 0.5) +  # x 방향 에러바
  labs(x = "Mean of STR length by LRS", y = "Mean of STR length by SRS") +
  theme_minimal() + 
  coord_fixed(ratio = 1) 

## 굿
eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% left_join(simple_concordance) %>% #head()
  filter(concordance_rate != 1) %>% 
  select(-concordance_point_sum,-concordance_rate) %>% #head(1000) %>%head()
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% #head(1000) %>% #head()
  ggplot(aes(x= mean_Whole.length_TRGT,y=mean_Whole.length_EH)) + 
  geom_hex(bins = 20) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "Mean of STR length by LRS", y = "Mean of STR length by SRS", fill = "# of STRs") +  # 축 및 범례 라벨
  coord_fixed(ratio = 1) +
  theme_step1() 


eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% left_join(simple_concordance) %>% #head()
  filter(concordance_rate != 1) %>% 
  select(-concordance_point_sum,-concordance_rate) %>% #head(1000) %>%head()
  pivot_wider(names_from = name,values_from = mean_Whole.length:min_Whole.length) %>% 
  select(STR_ID:mean_Whole.length_TRGT) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.notconrconrdacne1.v2.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  
  

library(ggplot2)
## 별로
eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% left_join(simple_concordance) %>% #head()
  filter(concordance_rate != 1) %>%  #head(1000) %>% #head()
  filter(name == "TRGT") %>%
  ggplot(aes(x=concordance_rate,y=mean_Whole.length)) + 
  geom_hex(bins = 15) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "Concordance rate", y = "mean of STR length", fill = "# of STRs") +  # 축 및 범례 라벨
  theme_step1()

  
eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength %>% left_join(simple_concordance) %>% #head()
  filter(concordance_rate != 1) %>%  #head(1000) %>% #head()
  filter(name == "TRGT") %>% head()
  

