library(tidyverse)

ref <- read.table("~/Desktop/KU/@research/STR/db/annovar/output_bed_ref")
ref %>% unique()-> ref

colnames(ref) <- c("chrom","start","end","STR_ID")

anno <- read_table("~/Desktop/KU/@research/STR/db/annovar/Raw.anno.merge.processing.onlyneed.txt")
anno %>% 
  mutate(type = case_when(
    type == "upstream;downstream" ~ "upstream",
    type == "UTR5;UTR3" ~ "UTR5",TRUE ~ type)) -> anno
concordance_rate <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt")

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR
head(common_STR)
common_STR %>% count(STR_DB %in% concordance_rate$STR_ID)
common_STR %>% filter(!(STR_DB %in% concordance_rate$STR_ID))

concordance_rate %>% filter(!(STR_ID %in% common_STR$STR_DB))

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro
head(final_ref_pro)
head(concordance_rate)
head(anno)

concordance_rate %>% left_join(anno) %>% left_join(ref) %>% left_join(final_ref_pro) -> df
head(df)

df %>% filter(concordance_rate != 1) %>% #count(type)
  group_by(RU.length,MOTIFS,GC,main_cyto_type,type) %>%  
  summarise(concordance_rate = mean(concordance_rate)) -> data
head(df)
###
eh_trgt_merge_simple_pass_intersect_forConcordance <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_trgt_merge_simple_pass_intersect_forConcordance.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance)

eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID:TRGT_STR2,RU:EH_STR2) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(concordance_point = STR1 + STR2) %>% group_by(STR_ID) %>% 
  summarise(concordance_point_sum = sum(concordance_point)) %>% 
  mutate(concordance_rate = concordance_point_sum/130) -> concordance_bySTR_simpleSTR

head(concordance_bySTR_simpleSTR)

#### without concordance == 1
eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID:TRGT_STR2,RU:EH_STR2) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(STR1_diff = TRGT_STR1 - EH_STR1,STR2_diff = TRGT_STR2 - EH_STR2) %>% #head()
  mutate(STR1_check = ifelse(STR1_diff < 0,"SRS",ifelse(STR1_diff == 0,"E","LRS")),STR2_check = ifelse(STR2_diff < 0,"SRS",ifelse(STR2_diff == 0,"E","LRS"))) %>%
  filter(STR_ID %in% concordance_rate[concordance_rate$concordance_rate != 1,]$STR_ID) -> concordance_bySTR_simpleSTR_v2


  
head(concordance_bySTR_simpleSTR_v2)
concordance_bySTR_simpleSTR_v2 %>% pivot_longer(STR1_check:STR2_check) %>% count(value)
#1 E     15291860
#2 LRS    1387889
#3 SRS    1693151
concordance_bySTR_simpleSTR_v2 %>% select(ID,STR_ID,TRGT_STR1,EH_STR1,STR1_diff,STR1_check) %>% mutate(allele = 1) -> concordance_bySTR_simpleSTR_v2_STR1
concordance_bySTR_simpleSTR_v2 %>% select(ID,STR_ID,TRGT_STR2,EH_STR2,STR2_diff,STR2_check) %>% mutate(allele = 2) -> concordance_bySTR_simpleSTR_v2_STR2
colnames(concordance_bySTR_simpleSTR_v2_STR1) <-c("ID","STR_ID","TRGT_STR","EH_STR","STR_diff","STR_check","allele")
colnames(concordance_bySTR_simpleSTR_v2_STR2) <-c("ID","STR_ID","TRGT_STR","EH_STR","STR_diff","STR_check","allele")
head(final_ref_pro)

concordance_bySTR_simpleSTR_v2_STR1 %>% rbind(concordance_bySTR_simpleSTR_v2_STR2) %>% 
  left_join(final_ref_pro) %>%
  mutate(STR_length = ifelse(STR_check == "SRS",EH_STR*RU.length,TRGT_STR*RU.length)) -> concordance_bySTR_simpleSTR_v2_STR_merge

head(concordance_bySTR_simpleSTR_v2_STR_merge)
concordance_bySTR_simpleSTR_v2_STR_merge %>% filter(STR_check == "SRS") %>% filter(STR_diff < 0)

concordance_bySTR_simpleSTR_v2_STR_merge %>% filter(STR_length > 100) %>% 
  ggplot(aes(x=RU.length,y=STR_length,color=STR_check)) +
  geom_point()

head(concordance_bySTR_simpleSTR_v2_STR_merge)
#write.table(concordance_bySTR_simpleSTR_v2_STR_merge,"~/Desktop/KU/@research/STR/figure/EH_TRGT_STR_length_compare_byID_merge.txt",col.names=T,row.names=F,quote=F,sep="\t")

###### all concordance


eh_trgt_merge_simple_pass_intersect_forConcordance %>% select(ID:TRGT_STR2,RU:EH_STR2) %>%
  mutate(STR1 = TRGT_STR1-EH_STR1,STR2 = TRGT_STR2-EH_STR2) %>% 
  mutate(STR1 = ifelse(STR1 == 0,1,0),STR2 = ifelse(STR2 == 0,1,0)) %>% #head()
  mutate(STR1_diff = TRGT_STR1 - EH_STR1,STR2_diff = TRGT_STR2 - EH_STR2) %>% #head()
  mutate(STR1_check = ifelse(STR1_diff < 0,"SRS",ifelse(STR1_diff == 0,"E","LRS")),STR2_check = ifelse(STR2_diff < 0,"SRS",ifelse(STR2_diff == 0,"E","LRS"))) -> 
  concordance_bySTR_simpleSTR_v2_allSTR




head(concordance_bySTR_simpleSTR_v2_allSTR)
concordance_bySTR_simpleSTR_v2_allSTR %>% pivot_longer(STR1_check:STR2_check) %>% count(value)
#1 E     37418250
#2 LRS    1387889
#3 SRS    1693151
concordance_bySTR_simpleSTR_v2_allSTR %>% select(ID,STR_ID,TRGT_STR1,EH_STR1,STR1_diff,STR1_check) %>% mutate(allele = 1) -> concordance_bySTR_simpleSTR_v2_allSTR_STR1
concordance_bySTR_simpleSTR_v2_allSTR %>% select(ID,STR_ID,TRGT_STR2,EH_STR2,STR2_diff,STR2_check) %>% mutate(allele = 2) -> concordance_bySTR_simpleSTR_v2_allSTR_STR2
colnames(concordance_bySTR_simpleSTR_v2_allSTR_STR1) <-c("ID","STR_ID","TRGT_STR","EH_STR","STR_diff","STR_check","allele")
colnames(concordance_bySTR_simpleSTR_v2_allSTR_STR2) <-c("ID","STR_ID","TRGT_STR","EH_STR","STR_diff","STR_check","allele")
head(final_ref_pro)

concordance_bySTR_simpleSTR_v2_allSTR_STR1 %>% rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR2) %>% 
  left_join(final_ref_pro) %>%
  mutate(STR_length = ifelse(STR_check == "SRS",EH_STR*RU.length,TRGT_STR*RU.length)) -> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge

head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% filter(STR_check == "SRS") %>% filter(STR_diff < 0)

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% filter(STR_length > 100) %>% 
  ggplot(aes(x=RU.length,y=STR_length,color=STR_check)) +
  geom_point()

head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)

## 이거는 Concordance == 1 allele 제거유무 파악
#write.table(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge,"~/Desktop/KU/@research/STR/figure/EH_TRGT_STR_length_compare_byID_merge_allSTR_withconcordance1.txt",col.names=T,row.names=F,quote=F,sep="\t")

#################
## 이거는 Concordance == 1 allele 제거유무 파악
concordance_bySTR_simpleSTR_v2_STR_merge <- read_table("~/Desktop/KU/@research/STR/figure/EH_TRGT_STR_length_compare_byID_merge.txt")
#concordance_bySTR_simpleSTR_v2_allSTR_STR_merge <- read_table("~/Desktop/KU/@research/STR/figure/EH_TRGT_STR_length_compare_byID_merge.txt")

head(concordance_bySTR_simpleSTR_v2_STR_merge)

## 이거는 Concordance == 1 allele 포함
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge <- read_table("~/Desktop/KU/@research/STR/figure/EH_TRGT_STR_length_compare_byID_merge_allSTR_withconcordance1.txt")
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)


#### 길이에 따른 concordance rate (long-read 길이 기준)
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% count(STR_length)
'''
STR_length        n
<chr>         <int>
  1 [0~50)     39215348
2 [100~150)     31083
3 [150~Inf)     42271
4 [50~100)    1210588
'''
#rm(concordance_model)
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(STR_length = "Overall") -> concordance_byID_allSTR

head(concordance_byID_allSTR)
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% group_by(STR_length) %>% count(STR_diff == 0) -> a
head(a)
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>% mutate(STR_length="Overall") %>%
  group_by(STR_length) %>% count(STR_diff == 0) -> a_overall
a_overall
  
a %>% rbind(a_overall) %>% rename(match = `STR_diff == 0`) %>% write.table("~/Desktop/KU/@research/STR/figure/figure2/f2.STR.allele.count.match.bySTR_length.txt",col.names = T,row.names = F,quote = F,sep = "\t")
STR.allele.count.match.bySTR_length <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.STR.allele.count.match.bySTR_length.txt")
STR.allele.count.match.bySTR_length %>% ggplot(aes(x=STR_length,y=n,fill=match)) + 
  geom_bar(stat = 'identity',position = "fill")
#a %>% pivot_wider(names_from = `STR_diff == 0`,values_from = n) -> a
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% group_by(ID,STR_length) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) -> a
head(a)

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% group_by(ID,STR_length) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% 
  select(ID,prop,STR_length) %>% rbind(concordance_byID_allSTR) %>% write.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance.range.byID_simpleSTR.v2.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")),
concordacen_byID <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance.range.byID_simpleSTR.v2.txt") 
concordacen_byID %>% ggplot(aes(x=factor(STR_length,levels=c("Overall","[0~50)","[50~100)","[100~150)","[150~Inf)")),
                                y=prop,fill=factor(STR_length,levels=c("Overall","[0~50)","[50~100)","[100~150)","[150~Inf)")))) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5)  + 
  labs(x="",y="") + 
  theme(legend.position = "none",
    legend.title = element_text()) + 
  theme_step1()
  
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(STR_length >= 50) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "upper50") -> concordance_byID_allSTR_upper50
head(concordance_byID_allSTR_upper50)


concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(STR_length >= 100) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "upper100") -> concordance_byID_allSTR_upper100

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(STR_length >= 150) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "upper150") -> concordance_byID_allSTR_upper150

concordance_byID_allSTR %>% rbind(concordance_byID_allSTR_upper50) %>% rbind(concordance_byID_allSTR_upper100) %>% 
  rbind(concordance_byID_allSTR_upper150) %>% ggplot(aes(x=g,y=prop,fill=g)) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  ylim(c(0,1))

## 길이 범위
#rm(eh_trgt_merge_simple_pass_intersect_forConcordance)

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(STR_length < 50) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "under50") -> concordance_byID_allSTR_under50
head(concordance_byID_allSTR_upper50)


concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(STR_length >= 50,STR_length<100) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "50~100") -> concordance_byID_allSTR_under50_100

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(STR_length >= 100,STR_length<150) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "100~150") -> concordance_byID_allSTR_under100_150

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(STR_length >= 150) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "150~") -> concordance_byID_allSTR_upper150


concordance_byID_allSTR_under50 %>% rbind(concordance_byID_allSTR_under50_100) %>% rbind(concordance_byID_allSTR_under100_150) %>% 
  rbind(concordance_byID_allSTR_upper150) %>% ggplot(aes(x=g,y=prop,fill=g)) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  ylim(c(0,1))


###
#### 길이에 따른 concordance rate repeat count 기준
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% 
  mutate(TRGT_STR = case_when(
    TRGT_STR < 5 ~ "[0~5)",
    TRGT_STR < 10 ~ "[5~10)",
    TRGT_STR < 20 ~ "[10~20)",
    TRUE ~ '[20~Inf)')) %>% count(TRGT_STR)

#1 [0~5)    12454620
#2 [10~20)   9796898
#3 [20~Inf)  2290020
#4 [5~10)   15957752

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "ALL") -> concordance_byID_allSTR

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% 
  filter(RU.length >= 5) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "upper5") -> concordance_byID_allSTR_RUupper5


concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(RU.length >= 10) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "upper10") -> concordance_byID_allSTR_RUupper10

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>%
  filter(RU.length >= 20) %>%  group_by(ID) %>% count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "upper20") -> concordance_byID_allSTR_RUupper20

concordance_byID_allSTR %>% rbind(concordance_byID_allSTR_RUupper5) %>% rbind(concordance_byID_allSTR_RUupper10) %>% 
  rbind(concordance_byID_allSTR_RUupper20) %>% ggplot(aes(x=g,y=prop,fill=g)) +
  geom_violin() + 
#  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  ylim(c(0,1))

### repeat count 범위
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% group_by(ID) %>% filter(RU.length < 5) %>%
  count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "<5") -> concordance_byID_allSTR_RU_0_5

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% group_by(ID) %>% filter(RU.length >= 5,RU.length < 10) %>%
  count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "5~10") -> concordance_byID_allSTR_RU_5_10
  

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% group_by(ID) %>% filter(RU.length >= 10,RU.length < 20) %>%
  count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "10~20") -> concordance_byID_allSTR_RU_10_20

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% group_by(ID) %>% filter(RU.length >= 20) %>%
  count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "20~") -> concordance_byID_allSTR_RU_20


concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% group_by(ID) %>% filter(RU.length >= 20) %>%
  count(STR_diff == 0) %>% 
  mutate(prop = prop.table(n)) %>% filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop) %>% mutate(g = "20~") -> concordance_byID_allSTR_RU_20



concordance_byID_allSTR_RU_0_5 %>% rbind(concordance_byID_allSTR_RU_5_10) %>% rbind(concordance_byID_allSTR_RU_10_20) %>% 
  rbind(concordance_byID_allSTR_RU_20) %>% ggplot(aes(x=g,y=prop,fill=g)) +
  geom_violin()
  #  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  #ylim(c(0,1))
####### GC
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>%  
  mutate(GC = case_when(
    GC < 0.25 ~ "[0~0.25)",
    GC < 0.5 ~ "[0.25~0.5)",
    GC < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(ID,GC) %>% count(STR_diff == 0) %>% #head()
  mutate(prop=prop.table(n)) %>%filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop,GC) %>% mutate(g = "1.ALL") -> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC
#memory.limit(size = 16000)  
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>% filter(STR_length < 50) %>%
  mutate(GC = case_when(GC < 0.25 ~ "[0~0.25)",GC < 0.5 ~ "[0.25~0.5)",GC < 0.75 ~ "[0.5~0.75)",TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(ID,GC) %>% count(STR_diff == 0) %>% #head()
  mutate(prop=prop.table(n)) %>%filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop,GC) %>% mutate(g = "2.<50")-> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_50

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>% filter(STR_length >= 50,STR_length < 100) %>%
  mutate(GC = case_when(GC < 0.25 ~ "[0~0.25)",GC < 0.5 ~ "[0.25~0.5)",GC < 0.75 ~ "[0.5~0.75)",TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(ID,GC) %>% count(STR_diff == 0) %>% #head()
  mutate(prop=prop.table(n)) %>%filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop,GC)%>%mutate(g = "3.<100") -> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_50_100

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>% filter(STR_length >= 100,STR_length < 150) %>%
  mutate(GC = case_when(GC < 0.25 ~ "[0~0.25)",GC < 0.5 ~ "[0.25~0.5)",GC < 0.75 ~ "[0.5~0.75)",TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(ID,GC) %>% count(STR_diff == 0) %>% #head()
  mutate(prop=prop.table(n)) %>%filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop,GC) %>% mutate(g = "4.<150")-> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_100_150

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = TRGT_STR * RU.length) %>% filter(STR_length >= 150) %>%
  mutate(GC = case_when(GC < 0.25 ~ "[0~0.25)",GC < 0.5 ~ "[0.25~0.5)",GC < 0.75 ~ "[0.5~0.75)",TRUE ~ '[0.75~1]')) %>% #count(GC)
  group_by(ID,GC) %>% count(STR_diff == 0) %>% #head()
  mutate(prop=prop.table(n)) %>%filter(`STR_diff == 0` == "TRUE") %>% select(ID,prop,GC) %>%mutate(g = "5.>=150") -> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_150


concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC %>% rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_50) %>%
  rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_50_100) %>% rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_100_150) %>%
  rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_150) %>% 
  mutate(STR_length = case_when(
    g == "1.ALL"  ~ "Overall",
    g == "2.<50"  ~ "[0~50)",
    g == "3.<100"  ~ "[50~100)",
    g == "4.<150" ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>% select(-g) %>% #head()
  write.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance.range.byID_length_GC_simpleSTR.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")),

concordance.range.byID_length_GC_simpleSTR <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance.range.byID_length_GC_simpleSTR.txt")
head(concordance.range.byID_length_GC_simpleSTR)
concordance.range.byID_length_GC_simpleSTR %>%#head()
  ggplot(aes(x=factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")),
           y=prop,fill=factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")))) + 
  geom_violin() + 
  labs(x="GC content (%)",y="Concordance") + 
  theme(legend.position = 'none') + 
  facet_grid(~STR_length)

concordance.range.byID_length_GC_simpleSTR.txt %>% 
  ggplot(aes(x=factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")),
             y=prop,fill=factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")))) + 
  geom_violin() + 
  labs(x="GC content (%) : ALL STR",y="Concordance") + 
  theme(legend.position = 'none')

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC %>% rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_50) %>%
  rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_50_100) %>% rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_100_150) %>%
  rbind(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_GC_150) %>% 
  ggplot(aes(x=factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")),
             y=prop,fill=factor(GC,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")))) + 
  geom_violin() + 
  labs(x="GC content (%)",y="Concordance") + 
  theme(legend.position = 'none') + 
  facet_grid(~g)
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_diff_logi = ifelse(TRGT_STR == EH_STR,1,0)) %>%
  mutate(STR_length = TRGT_STR * RU.length) -> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% mutate(STR_length = rescale(STR_length)) %>% select(STR_length,STR_diff_logi,GC,RU.length,TRGT_STR)-> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge
rm(a)
library(scales)

#glm(STR_diff_logi ~ STR_length + GC + RU.length,data=concordance_bySTR_simpleSTR_v2_allSTR_STR_merge,family = binomial) -> concordance_model
glm(STR_diff_logi ~ STR_length,data=concordance_bySTR_simpleSTR_v2_allSTR_STR_merge,family = binomial) -> concordance_model
glm(STR_diff_logi ~ GC,data=concordance_bySTR_simpleSTR_v2_allSTR_STR_merge,family = binomial) -> concordance_model_GC
glm(STR_diff_logi ~ RU.length,data=concordance_bySTR_simpleSTR_v2_allSTR_STR_merge,family = binomial) -> concordance_model_RU.length
glm(STR_diff_logi ~ TRGT_STR,data=concordance_bySTR_simpleSTR_v2_allSTR_STR_merge,family = binomial) -> concordance_model_MC
glm(STR_diff_logi ~ STR_length + GC + RU.length + TRGT_STR,data=concordance_bySTR_simpleSTR_v2_allSTR_STR_merge,family = binomial) -> concordance_model_all
summary(concordance_model)

# 계수 출력
coefficients(concordance_model)

# 계수의 지수값 (오즈비로 변환)
exp(coefficients(concordance_model))

#•	Intercept (절편): exp(3.459727) = 31.80829
#•	STR_length: exp(-474.423539) = 9.130084e-207 (즉, 거의 0에 가까움)

summary(concordance_model_MC)
'''
Coefficients:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.700e+00  9.904e-04  1716.9   <2e-16 ***
  TRGT_STR    -8.344e-03  6.482e-05  -128.7   <2e-16 ***
'''
head(STR_diff_logi)
coefficients(concordance_model_GC)
exp(coefficients(concordance_model_GC))

table(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge$STR_diff_logi)

'''
(Intercept)  1.422e+00  6.640e-04  2141.9   <2e-16 ***
  STR_length  -2.355e+02  1.200e-01 -1962.0   <2e-16 ***
  GC           2.214e+00  8.716e-04  2540.2   <2e-16 ***
  RU.length    7.264e-02  1.965e-04   369.7   <2e-16 ***
  '''
coefficients(concordance_model_all)
exp(coefficients(concordance_model_all))


'''
  GC                n
<chr>         <int>
  1 [0.25~0.5)  8460010
2 [0.5~0.75) 15375880
3 [0.75~1]    1680900
4 [0~0.25)   14982500
'''


concordance_byID_allSTR_RU_0_5 %>% rbind(concordance_byID_allSTR_RU_5_10) %>% rbind(concordance_byID_allSTR_RU_10_20) %>% 
  rbind(concordance_byID_allSTR_RU_20) %>% ggplot(aes(x=g,y=prop,fill=g)) +
  geom_violin()


memory.limit()





################
## 확인
concordance_bySTR_simpleSTR_v2_STR_merge %>% mutate(STR_length_diff = STR_diff*RU.length) %>% filter(STR_check != "E") %>%
  filter(abs(STR_length) >= 100) %>%
#  filter(abs(STR_length_diff) >= 100) %>%
  ggplot(aes(x=RU.length,y=STR_length_diff,color=STR_check)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1) 

## 확인 -> 이건 의미가 있을듯 by STR length STR_length
head(concordance_bySTR_simpleSTR_v2_STR_merge)
concordance_bySTR_simpleSTR_v2_STR_merge %>% mutate(STR_length_diff = STR_diff*RU.length) %>% #filter(STR_check != "E") %>%
  filter(abs(STR_length) >= 100) %>% #head()
  ggplot(aes(x=STR_length,y=STR_length_diff,color=STR_check)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1) #+
#  theme_step1()

## G 는 별로 
head(concordance_bySTR_simpleSTR_v2_STR_merge)
concordance_bySTR_simpleSTR_v2_STR_merge %>% mutate(STR_length_diff = STR_diff*RU.length) %>% #filter(STR_check != "E") %>%
  filter(abs(STR_length) >= 100) %>% #head()
  ggplot(aes(x=GC,y=STR_length_diff,color=STR_check)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1)
#  theme_step1()
#bin_genome()
concordance_rate %>% left_join(final_ref_pro) %>% #filter(concordance_rate != 1) %>%
  ggplot(aes(x=GC,y=concordance_rate)) + 
  geom_hex(bins=15) + 
  geom_smooth(method = "lm", se = TRUE)

concordance_rate %>% left_join(final_ref_pro) %>% #filter(concordance_rate != 1) %>%
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
  concordance_rate == 1 ~ "1")) %>% 
  ggplot(aes(x=concordance_rate_range,y=GC,fill=concordance_rate_range)) + 
  geom_boxplot()

concordance_rate %>% left_join(final_ref_pro) %>% #head()
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
    concordance_rate == 1 ~ "1")) %>% group_by(concordance_rate_range) %>% summarise(mean_GC = mean(GC)) %>%
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                       "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                       "[0.8,0.9)", "[0.9,1)","1")),y=mean_GC)) + 
  geom_point() + 
  labs(x="Concordance Interval",y="Mean of GC contents")



## coverage 
#eh_simple_LC <- read_table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_LC.txt") %>% filter(STR_ID %in% common_STR$STR_DB)
eh_simple_LC <- read_table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_oribam_coverage.txt") %>% filter(STR_ID %in% common_STR$STR_DB)
trgt_simple_coverage <- read_table("~/Desktop/KU/@research/STR/trgt/trgt_simpleSTR_coverage.txt") #%>% select(ID,STR_ID,meandepth)
trgt_simple_coverage %>% select(ID,STR_ID,meandepth)
head(eh_simple_LC)
head(trgt_simple_coverage)
trgt_simple_coverage %>% filter(!str_detect(STR_ID,"chr"))
head(concordance_bySTR_simpleSTR_v2_STR_merge)

## 별 차이 없음
concordance_bySTR_simpleSTR_v2_STR_merge %>% #filter(STR_check != "E") %>% 
  left_join(trgt_simple_coverage) %>% #head()
  ggplot(aes(x=STR_check,y=meandepth,fill=STR_check)) + 
  geom_violin()

head(concordance_rate)
head(eh_simple_LC)
eh_simple_LC %>% filter(LC > 50000)
eh_simple_LC %>% select(STR_ID) %>% unique %>% dim()
eh_simple_LC %>% select(STR_ID) %>% filter(STR_ID == "ZNF713")
trgt_simple_coverage %>% select(STR_ID) %>% filter(STR_ID == "ZNF713")
eh_simple_LC %>% select(STR_ID) %>% filter(!(STR_ID %in% concordance_rate$STR_ID))

trgt_simple_coverage %>% select(STR_ID) %>% unique %>% dim()
trgt_simple_coverage %>% filter(!(STR_ID %in% concordance_rate$STR_ID))
eh_simple_LC %>% left_join(trgt_simple_coverage) %>% left_join(concordance_rate) %>% 
  mutate(g = ifelse(concordance_rate == 1,"1","0")) %>%
  group_by(STR_ID,g) %>% head()
  summarise(mean_LC = mean(LC),mean_trgt = mean(meandepth)) -> concordance_depth

head(concordance_depth)
trgt_simple_coverage %>% filter(STR_ID == "ZNF713")
#concordance_depth %>% filter(is.na(mean_trgt))

ggplot(concordance_depth,aes(x=mean_LC,y=mean_trgt)) + 
  geom_hex(bins=30) +
  facet_grid(~g)
 

concordance_depth %>% left_join(concordance_rate) %>% #filter(concordance_rate == 1) %>%
  ggplot(aes(x=mean_trgt,y=concordance_rate)) + 
  geom_point() +
  labs(x="Mean of mapping depth by samtools coverage (LRS)") + 
  geom_smooth(method = "lm", se = TRUE)


concordance_depth %>% left_join(concordance_rate) %>% #filter(concordance_rate == 1) %>%
  ggplot(aes(x=mean_LC,y=concordance_rate)) + 
  geom_point() +
  labs(x="Mean of Locus Coverage by EH VCF (SRS)") + 
  geom_smooth(method = "lm", se = TRUE)




concordance_depth %>% left_join(concordance_rate) %>% filter(concordance_rate == 0) %>%
  ggplot(aes(x=mean_trgt,y=mean_LC)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)



head(trgt_simple_coverage)

concordance_bySTR_simpleSTR_v2_STR_merge %>% mutate(STR_length_diff = STR_diff*RU.length) %>% filter(STR_check != "E") %>%
  #filter(abs(STR_length) >= 100) %>%
  filter(abs(STR_length_diff) >= 10000) %>% head()

head(concordance_bySTR_simpleSTR)
head(concordance_bySTR_simpleSTR_v2_STR_merge)
concordance_bySTR_simpleSTR_v2_STR_merge %>% filter(STR_length >= 100) %>% left_join(concordance_bySTR_simpleSTR) %>% count(STR)

##의미 없음
concordance_bySTR_simpleSTR_v2_STR_merge %>% #head()
  mutate(STR_length_range = case_when(
    STR_length < 10 ~ "[0,10)",
    STR_length >= 10 & STR_length < 50 ~ "[10,50)",
    STR_length >= 50 & STR_length < 100 ~ "[50,100)",
    STR_length >= 100 & STR_length < 200 ~ "[100,200)",
    STR_length >= 200 & STR_length < 300 ~ "[200,300)",
    STR_length >= 300 & STR_length < 400 ~ "[300,400)",
    STR_length >= 400 & STR_length < 500 ~ "[400,500)",
    STR_length >= 500 & STR_length < 1000 ~ "[500,1000)",
    TRUE ~ "1000+")) %>% #head()
  ggplot(aes(x=STR_length_range,y=STR_diff)) + 
  geom_violin()
  
library(ggplot2)

## samtools ori bam shortread
eh_simple_coverage <- read_table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_oribam_coverage.txt") %>% filter(STR_ID %in% common_STR$STR_DB)
trgt_simple_coverage <- read_table("~/Desktop/KU/@research/STR/trgt/trgt_simpleSTR_coverage.txt") #%>% select(ID,STR_ID,meandepth)

head(eh_simple_coverage)
head(trgt_simple_coverage)
head(concordance_rate)

eh_simple_coverage %>% group_by(STR_ID) %>% summarise(meandepth = mean(meandepth),meanbaseq = mean(meanbaseq),meanmapq = mean(meanmapq)) -> eh_simple_coverage_mean
trgt_simple_coverage %>% group_by(STR_ID) %>% summarise(meandepth = mean(meandepth),meanbaseq = mean(meanbaseq),meanmapq = mean(meanmapq)) -> trgt_simple_coverage_mean

head(eh_simple_coverage_mean)
head(trgt_simple_coverage_mean)
colnames(eh_simple_coverage_mean) <- c("STR_ID","EH_meandepth","EH_meanbaseq","EH_maenmapq")
colnames(trgt_simple_coverage_mean) <- c("STR_ID","trgt_meandepth","trgt_meanbaseq","trgt_maenmapq")
trgt_simple_coverage_mean

concordance_rate %>% left_join(trgt_simple_coverage_mean) %>% left_join(eh_simple_coverage_mean) %>%
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
    concordance_rate == 1 ~ "1")) -> merge_coverage

merge_coverage %>% 
  ggplot(aes(x=trgt_meandepth,y=EH_meandepth)) +
  geom_point() + 
  facet_grid(~concordance_rate_range)
  
merge_coverage %>% 
  ggplot(aes(x=trgt_meanbaseq,y=EH_meanbaseq)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # 회귀선 추가 (선형 회귀)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  ylim(c(0,60)) + xlim(c(0,60)) + 
  labs(x="mean baseQ (LRS)",y="mean baseQ (SRS)") + 
  facet_grid(~factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                     "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                     "[0.8,0.9)", "[0.9,1)","1")))

merge_coverage %>% 
  ggplot(aes(x=trgt_maenmapq,y=EH_maenmapq)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # 회귀선 추가 (선형 회귀)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  ylim(c(0,60)) + 
  labs(x="mean mapQ (LRS)",y="mean mapQ (SRS)") + 
  facet_grid(~factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                      "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                      "[0.8,0.9)", "[0.9,1)","1")))


merge_coverage %>% filter(concordance_rate == 0) %>%
  ggplot(aes(x=EH_maenmapq,y=EH_maenmapq)) +
  geom_point() + 
  labs(x="mean mapQ (LRS)",y="mean mapQ (SRS)") 


head(merge_coverage)
str_length <- read.table()
merge_coverage %>% 
  mutate(concordance_rate_range = case_when(
    concordance_rate < 0.4 ~ "[0, 0.4)",
    concordance_rate < 0.8 ~ "[0.4, 0.8)",
    TRUE ~ "[0.8, 1]"
  )) %>%
  ggplot(aes(x=trgt_maenmapq,y=EH_maenmapq,alpha=0.5)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # 회귀선 추가 (선형 회귀)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  ylim(c(0,60)) + 
  labs(x="mean mapQ (LRS)",y="mean mapQ (SRS)") + 
  facet_grid(~factor(concordance_rate_range,levels= c("[0, 0.4)","[0.4, 0.8)","[0.8, 1]")))


#write.table(merge_coverage,"~/Desktop/KU/@research/STR/figure/extra_info/STR.mean.bam.stats.txt",col.names = T,row.names = F,quote = F,sep = "\t")


################## methylation
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge <- read_table("~/Desktop/KU/@research/STR/figure/EH_TRGT_STR_length_compare_byID_merge_allSTR_withconcordance1.txt")

trgt_simple_pass_intersect_merge_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure2/trgt_simple_pass_intersect_merge_processing.txt")
head(trgt_simple_pass_intersect_merge_processing)
trgt_simple_pass_intersect_merge_processing %>% select(TRID) %>% unique() %>% dim()

trgt_simple_pass_intersect_merge_processing %>% filter(TRGT_AM1 != -1 & TRGT_AM2 != -1) %>% count(TRID) %>% filter(n == 65) -> trgt_simple_pass_intersect_merge_processing_allMethyl

trgt_simple_pass_intersect_merge_processing %>% 
  filter(TRGT_AM1 != -1 & TRGT_AM2 != -1) %>% count(TRID) %>% filter(n == 65) %>% dim()

trgt_simple_pass_intersect_merge_processing %>% 
  filter(TRGT_AM1 != -1 & TRGT_AM2 != -1) %>% count(TRID) %>% filter(n >= 59) %>% dim()


head(trgt_simple_pass_intersect_merge_processing_allMethyl)
trgt_simple_pass_intersect_merge_processing_allMethyl %>% count(TRID %in%common_STR$STR_DB)
head(trgt_simple_pass_intersect_merge_processing)
head(trgt_simple_pass_intersect_merge_processing_allMethyl)
trgt_simple_pass_intersect_merge_processing_allMethyl %>% select(TRID) %>% dim()

trgt_simple_pass_intersect_merge_processing %>% select(ID,TRID,TRGT_AM1) %>% mutate(allele = 1) %>% filter(TRID %in% trgt_simple_pass_intersect_merge_processing_allMethyl$TRID) -> trgt_simple_pass_intersect_merge_processing_STR1
trgt_simple_pass_intersect_merge_processing %>% select(ID,TRID,TRGT_AM2) %>% mutate(allele = 2) %>% filter(TRID %in% trgt_simple_pass_intersect_merge_processing_allMethyl$TRID) -> trgt_simple_pass_intersect_merge_processing_STR2

colnames(trgt_simple_pass_intersect_merge_processing_STR1) <-c("ID","STR_ID","TRGT_AM","allele")
colnames(trgt_simple_pass_intersect_merge_processing_STR2) <-c("ID","STR_ID","TRGT_AM","allele")

head(trgt_simple_pass_intersect_merge_processing_allMethyl)
trgt_simple_pass_intersect_merge_processing_STR1 %>% rbind(trgt_simple_pass_intersect_merge_processing_STR2) %>% filter(STR_ID %in% common_STR$STR_DB) -> trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall

head(trgt_simple_pass_intersect_merge_processing_STR)
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)
#concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% filter(STR_ID %in% trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall$STR_ID)
head(trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall)
head(trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall)

trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>% count(STR_ID %in% concordance_bySTR_simpleSTR_v2_allSTR_STR_merge$STR_ID)
trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>% left_join(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge) -> trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall

head(trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall)

#write.table(trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall,"~/Desktop/KU/@research/STR/figure/trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall.txt",col.names = T,row.names = F,quote = F,sep = "\t")
trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall <-read_table("~/Desktop/KU/@research/STR/figure/trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall.txt")

head(trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall)
dim(trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall)

## 별로
trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>%
  ggplot(aes(x=TRGT_AM,y=STR_diff,fill=STR_check)) + 
  geom_point()

head(trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall)
trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>% select(STR_ID) %>% unique() %>% dim()

trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>% left_join(concordance_rate) %>% 
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
    concordance_rate == 1 ~ "1")) %>% select(STR_ID,concordance_rate_range) %>% unique() %>% count(concordance_rate_range)

## 별로
trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>% left_join(concordance_rate) %>%
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
    concordance_rate == 1 ~ "1")) %>% group_by(STR_ID,concordance_rate_range) %>% summarise(mean_AM = mean(TRGT_AM)) %>%
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                       "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                       "[0.8,0.9)", "[0.9,1)","1")),y=mean_AM)) + 
  geom_violin()


## 별로
trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>% left_join(concordance_rate) %>% 
  ggplot(aes(x=TRGT_AM,y=concordance_rate)) + 
  geom_hex(bins=15) +
  geom_smooth(method = "lm", se = FALSE)


trgt_simple_pass_intersect_merge_processing_STR_onlyallMethylcall %>% left_join(concordance_rate) %>%
  write.table("~/Desktop/KU/@research/STR/figure/methylation.allsamplecallremain.STRQCinfo.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  


#################
concordance_bySTR_simpleSTR_v2 %>% select(ID,STR_ID,STR1_check,STR2_check) %>%
  pivot_longer(STR1_check:STR2_check) %>% count(STR_ID,value) -> concordance_bySTR_simpleSTR_v2_check

head(concordance_bySTR_simpleSTR_v2_check)
concordance_bySTR_simpleSTR_v2_check %>% filter(n == 130) -> concordance_bySTR_simpleSTR_v2_check_allsample_overestiname

head(concordance_bySTR_simpleSTR_v2_check_allsample_overestiname)

trgt_simple_pass_intersect_merge_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure2/trgt_simple_pass_intersect_merge_processing.txt")
head(trgt_simple_pass_intersect_merge_processing)
eh_simple_pass_intersect_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_simple_pass_intersect_processing.txt")

ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)
ref %>% select(Revio,Illumina) -> ref
head(ref)
colnames(ref) <- c("ID","EH_ID")

head(concordance_bySTR_simpleSTR_v2_check_allsample_overestiname)
concordance_bySTR_simpleSTR_v2_check_allsample_overestiname %>% filter(value == "SRS") -> concordance_bySTR_simpleSTR_v2_check_allsample_overestiname_SRS
concordance_bySTR_simpleSTR_v2_check_allsample_overestiname %>% filter(value == "LRS") -> concordance_bySTR_simpleSTR_v2_check_allsample_overestiname_LRS

#eh_simple_pass_intersect_processing %>% count(ID) %>% dim() #[1] 65  2
head(eh_simple_pass_intersect_processing)
head(concordance_bySTR_simpleSTR_v2_STR_merge)
eh_simple_pass_intersect_processing  %>% rename(EH_ID = ID) %>% left_join(ref) %>% select(ID,STR_ID,EH_type1:EH_type2) %>% #head()
  filter(STR_ID %in% concordance_bySTR_simpleSTR_v2_check_allsample_overestiname_SRS$STR_ID) %>% pivot_longer(EH_type1:EH_type2) %>% #head()
  mutate(allele = ifelse(name == "EH_type1",1,2)) %>% #head()
  left_join(concordance_bySTR_simpleSTR_v2_STR_merge %>% select(ID,STR_ID,allele,EH_STR,STR_diff,GC)) -> eh_simple_pass_intersect_processing_type_check

head(eh_simple_pass_intersect_processing_type_check)
eh_simple_pass_intersect_processing_type_check %>% 
  ggplot(aes(x=EH_STR,y=abs(STR_diff),color=GC)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1) 

#FLANKING    52
#INREPEAT    10
#SPANNING 12548


head(concordance_bySTR_simpleSTR)



#####
library(ComplexHeatmap)
library(circlize)  # colorRamp2 사용을 위해 필요

methylation.allsamplecallremain.STRQCinfo.txt <- read_table("~/Desktop/KU/@research/STR/figure/methylation.allsamplecallremain.STRQCinfo.txt")
head(methylation.allsamplecallremain.STRQCinfo.txt)

anno <- read_table("~/Desktop/KU/@research/STR/db/annovar/Raw.anno.merge.processing.onlyneed.txt") %>% select(STR_ID,main_cyto_type,type)

head(anno)
methylation.allsamplecallremain.STRQCinfo.txt %>% arrange(ID,STR_ID,allele) %>%  
  mutate(STR_length_diff = STR_diff*RU.length) %>% filter(concordance_rate != 1) %>%
  left_join(anno)-> data
head(data)
dim(data)
# 필요한 패키지 로드
library(ComplexHeatmap)

head(data)
forHeat_toy<- data  %>% mutate(ID = paste0(ID,"_",allele))
head(forHeat_toy)
forHeat_toy %>% select(ID,STR_ID,TRGT_AM) %>% pivot_wider(names_from = STR_ID,values_from = TRGT_AM) -> forHeat_toy
head(forHeat_toy)
forHeat_toy_m <- as.matrix(forHeat_toy)
head(forHeat_toy_m)
forHeat_toy_m[1,1:5]
forHeat_toy_m <-apply(forHeat_toy_m, 2, as.numeric)
rownames(forHeat_toy_m) <- forHeat_toy$ID
forHeat_toy_m[,2:ncol(forHeat_toy_m)] -> forHeat_toy_m
colnames(forHeat_toy_m)
str_length(colnames(forHeat_toy_m))
#str_length(gsub("[^CG]", "",colnames(forHeat_toy_m)))/str_length(colnames(forHeat_toy_m))
head(forHeat_toy)
dim(forHeat_toy)
head(data)
data %>% select(STR_ID,GC) %>% unique() %>% dim()
data %>% select(STR_ID) %>% unique() %>% left_join(final_ref_pro) %>% select(STR_ID,RU.length) -> df_rulength
data %>% select(STR_ID) %>% unique() %>% left_join(final_ref_pro) %>% select(STR_ID,GC) -> df_GC
data %>% select(STR_ID,concordance_rate) %>% unique() -> df_concordance
#data %>% select(STR_ID,STR_diff) %>% group_by(STR_ID) %>% summarise(STR_diff = mean(abs(STR_diff))) %>%
 # mutate(norm_STR_diff = scale(STR_diff))-> df_STR_diff
data %>% mutate(STR_length = RU.length*TRGT_STR) %>% 
  group_by(STR_ID) %>% summarise(STR_length = mean(STR_length)) %>% #filter(STR_length < 50) %>% dim()
  mutate(STR_length_range = case_when(
    STR_length < 20 ~ "20",
    STR_length < 30 ~ "30",
    STR_length < 40 ~ "40",
    STR_length < 50 ~ "50",
    TRUE ~ "long"))-> df_str_length

library(RColorBrewer)
library(circlize)

#levels <- unique(df_str_length$STR_length_range)
levels <- c("long","50","40","30","20")
colors <- colorRampPalette(brewer.pal(5, "Greens"))(length(levels))
colors
# STR_length_range와 색상 매핑
col_map <- setNames(colors, rev(levels))


data %>% select(STR_ID,type,main_cyto_type) %>% unique()-> df_STR_anno

head(df_GC)
head(df_rulength)
head(df_concordance)
gpar(fill = 2:4)
anno_block(gp = gpar(fill = 2:4))
# Create annotations for RU.length and GC content
#col_anno <- columnAnnotation(concordance = df_concordance$concordance_rate,df_STR_diff = df_STR_diff$norm_STR_diff,RU.length = df_rulength$RU.length, GC = df_GC$GC)
col_anno = HeatmapAnnotation(
  concordance = df_concordance$concordance_rate,
  #df_STR_diff = df_STR_diff$norm_STR_diff,
  STR_length = df_str_length$STR_length_range, 
  GC = df_GC$GC,
  anno_type = df_STR_anno$type,
  anno_cyto = df_STR_anno$main_cyto_type,
  col = list(concordance = colorRamp2(c(0, 0.5,1), c("blue",'grey' ,"red")),
             GC = colorRamp2(c(0, 1), c("white","blue")),
             STR_length = col_map),
  simple_anno_size = unit(1, "cm"),
  annotation_name_side = "left"
)
## 되는것

# Generate the heatmap with clustering
Heatmap(forHeat_toy_m, 
        name = "methyl", 
        #col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_dend = TRUE, 
        show_column_dend = TRUE, 
        row_dend_reorder = TRUE,
        show_column_names = FALSE,  # Hide x-axis text (sample names)
        show_row_names = FALSE,  
        column_dend_reorder = TRUE,
        top_annotation = col_anno)



# 데이터 예시
# 필요한 패키지 로드
library(ComplexHeatmap)
library(circlize)
library(reshape2)
head(data)
anyNA(forHeat_toy)
forHeat_toy<- data  %>% mutate(ID = paste0(ID,"_",allele)) %>% filter(concordance_rate != 1)
forHeat_toy %>% select(ID,STR_ID,STR_length_diff) %>% pivot_wider(names_from = STR_ID,values_from = STR_length_diff) -> forHeat_toy
forHeat_toy_m <- as.matrix(forHeat_toy)
anyNA(forHeat_toy_m)

#forHeat_toy_m <-apply(forHeat_toy_m, 2, as.numeric)
rownames(forHeat_toy_m) <- forHeat_toy$ID
anyNA(forHeat_toy_m)

forHeat_toy_m[,2:ncol(forHeat_toy_m)] -> forHeat_toy_m
anyNA(forHeat_toy_m)
forHeat_toy_m <-apply(forHeat_toy_m, 2, as.double)
anyNA(forHeat_toy_m)

#forHeat_toy_m <- scale(forHeat_toy_m)
forHeat_toy_m[1:5,1:5]
anyNA(forHeat_toy_m)
forHeat_toy_m[is.na(forHeat_toy_m)] 
head(data)
data %>% select(STR_ID) %>% unique() %>% left_join(final_ref_pro) %>% select(STR_ID,GC) -> df_GC
data %>% select(STR_ID,concordance_rate) %>% unique() -> df_concordance
data %>% select(STR_ID,TRGT_AM) %>% group_by(STR_ID) %>% 
  summarise(mean_AM = mean(TRGT_AM)) -> df_am

library(RColorBrewer)
col_anno = HeatmapAnnotation(
  concordance = df_concordance$concordance_rate,
  GC = df_GC$GC,
  methyl = df_am$mean_AM,
  col = list(concordance = colorRamp2(c(0, 0.5,1), c("red",'grey' ,"blue")),
             GC = colorRamp2(c(0, 1), c("white","blue")),
             methyl = colorRamp2(c(0, 1), c("white","brown"))),
  simple_anno_size = unit(1, "cm"),
  annotation_name_side = "left"
)
## 되는것
head(data)
str(data$STR_length_diff)
summary(data$STR_length_diff)

data %>% summarise(min(STR_length_diff),max(STR_length_diff))
# Generate the heatmap with clustering
Heatmap(forHeat_toy_m, 
        name = "STR length dff",
        #col = colorRamp2(c(min(data$STR_length_diff),2 ,0, 2,max(data$STR_length_diff)), c("purple","blue", "white", "green","red")), 
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_dend = TRUE, 
        show_column_dend = TRUE, 
        row_dend_reorder = TRUE,
        show_column_names = FALSE,  # Hide x-axis text (sample names)
        show_row_names = FALSE,  
        column_dend_reorder = TRUE,
        top_annotation = col_anno)


# 히트맵 그리기
draw(ht_methylation)





head(data)
head(forHeat_toy)
data %>% select(STR_ID) %>% unique() %>% left_join(final_ref_pro) %>% select(STR_ID,GC) -> df_GC
data %>% select(STR_ID,concordance_rate) %>% unique() -> df_concordance
data %>% select(STR_ID,TRGT_AM) %>% group_by(STR_ID) %>% 
  summarise(mean_AM = mean(TRGT_AM)) -> df_am

library(RColorBrewer)
col_anno = HeatmapAnnotation(
  concordance = df_concordance$concordance_rate,
  GC = df_GC$GC,
  col = list(concordance = colorRamp2(c(0, 0.5,1), c("blue",'grey' ,"red")),
             GC = colorRamp2(c(0, 1), c("white","purple"))),
  simple_anno_size = unit(1, "cm"),
  annotation_name_side = "left"
)

####### methyl, STR_length, STR_diff_length
forHeat_toy<- data  %>% mutate(ID = paste0(ID,"_",allele))
forHeat_toy %>% select(ID,STR_ID,TRGT_AM) %>% pivot_wider(names_from = STR_ID,values_from = TRGT_AM) -> forHeat_toy
forHeat_toy_am <- as.matrix(forHeat_toy)
anyNA(forHeat_toy_am)
rownames(forHeat_toy_am) <- forHeat_toy$ID
anyNA(forHeat_toy_am)

forHeat_toy_am[,2:ncol(forHeat_toy_am)] -> forHeat_toy_am
anyNA(forHeat_toy_am)
forHeat_toy_am <-apply(forHeat_toy_am, 2, as.double)
anyNA(forHeat_toy_am)

#forHeat_toy_m <- scale(forHeat_toy_m)
forHeat_toy_am[1:5,1:5]
anyNA(forHeat_toy_am)
forHeat_toy_am[is.na(forHeat_toy_am)] 
head(data)



## 되는것
# Generate the heatmap with clustering
Heatmap(forHeat_toy_am, 
        name = "Methyl",
        #col = colorRamp2(c(min(data$STR_length_diff),2 ,0, 2,max(data$STR_length_diff)), c("purple","blue", "white", "green","red")), 
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_dend = TRUE, 
        show_column_dend = TRUE, 
        row_dend_reorder = TRUE,
        show_column_names = FALSE,  # Hide x-axis text (sample names)
        show_row_names = FALSE,  
        column_dend_reorder = TRUE,
        top_annotation = col_anno)

###
head(data)
forHeat_toy<- data  %>% mutate(ID = paste0(ID,"_",allele))
forHeat_toy %>% select(ID,STR_ID,STR_length) %>% pivot_wider(names_from = STR_ID,values_from = STR_length) -> forHeat_toy
forHeat_toy_strlength <- as.matrix(forHeat_toy)
anyNA(forHeat_toy_strlength)
rownames(forHeat_toy_strlength) <- forHeat_toy$ID
anyNA(forHeat_toy_strlength)

forHeat_toy_strlength[,2:ncol(forHeat_toy_strlength)] -> forHeat_toy_strlength
anyNA(forHeat_toy_strlength)
forHeat_toy_strlength <-apply(forHeat_toy_strlength, 2, as.double)
anyNA(forHeat_toy_strlength)

#forHeat_toy_m <- scale(forHeat_toy_m)
forHeat_toy_am[1:5,1:5]
anyNA(forHeat_toy_strlength)
forHeat_toy_am[is.na(forHeat_toy_strlength)] 
head(data)



## 되는것
# Generate the heatmap with clustering
Heatmap(forHeat_toy_strlength, 
        name = "STRlength",
        #col = colorRamp2(c(min(data$STR_length_diff),2 ,0, 2,max(data$STR_length_diff)), c("purple","blue", "white", "green","red")), 
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_dend = TRUE, 
        show_column_dend = TRUE, 
        row_dend_reorder = TRUE,
        show_column_names = FALSE,  # Hide x-axis text (sample names)
        show_row_names = FALSE,  
        column_dend_reorder = TRUE,
        top_annotation = col_anno) #-> a

row_clusters <- row_order(a)
print(row_clusters)

column_clusters <- column_order(a)
print(column_clusters)

row.names(forHeat_toy_strlength)[row_clusters]

head(data)
row_dend <- row_dend(a)
row_clusters <- cutree(as.hclust(row_dend), k = 2)  # 2개의 그룹으로 나누기
print(row_clusters) 
head(forHeat_toy)
forHeat_toy[1:5,1:5]
original_row_names <- rownames(forHeat_toy_strlength)
print(original_row_names)
row_clusters %>% as.data.frame()
data.frame(index=row_clusters,aID=forHeat_toy$ID) -> original_row_names
gender <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/revio.sex.info.txt")
data %>% select(ID) %>% left_join(gender) 

original_row_names %>% mutate(ID = str_split_fixed(aID,"_",2)[,1]) %>% 
  left_join(gender) %>% count(index,sex)

## homo diner
concordance_rate <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt")
head(concordance_rate)
head(final_ref)
library(tidyverse)
concordance_rate %>% left_join(final_ref %>% mutate(STR_ID = ID) %>% select(MOTIFS,STR_ID)) %>%
  mutate(MOTIF_type = case_when(
    grepl("^(A+|T+|C+|G+)$", MOTIFS) ~ "homopolymer",  # 단일 염기의 반복
    grepl("^(AT|TA|CG|GC|GT|TG|AG|GA|CT|TC|AC|CA)+$", MOTIFS) & nchar(MOTIFS) == 2 ~ "dinucleotide",  # 2개 염기 반복
    nchar(MOTIFS) == 3 ~ "trinucleotide",  # 3개 염기 반복
    nchar(MOTIFS) == 4 ~ "tetranucleotide",  # 4개 염기 반복
    TRUE ~ "other"  # 그 외의 경우  # 그 외의 경우
  )) %>% #filter(MOTIF_type == "other") %>% select(MOTIFS) %>% list()
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
    concordance_rate == 1 ~ "1")) %>% na.omit() %>% #count(MOTIF_type)
 count(MOTIF_type,concordance_rate_range) -> a
a %>%
  ggplot(aes(x=concordance_rate_range,y=n,fill=factor(MOTIF_type,levels=c("homopolymer",'dinucleotide',"trinucleotide","tetranucleotide","other")))) + 
  geom_bar(stat = "identity",position = "fill") + 
  theme_bw() + 
  theme(legend.title = element_blank())

a %>% filter(MOTIF_type == "homopolymer") %>%
  ggplot(aes(x=concordance_rate_range,y=n)) +
  labs(y="# of homopolymer STR") + 
  geom_bar(stat = "identity")



