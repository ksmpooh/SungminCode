### main figure 3 data processing # cirud plot and aano
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


eh_trgt_merge_simple_pass_intersect_forConcordance <- read_table("~/Desktop/KU/@research/STR/figure/figure2/eh_trgt_merge_simple_pass_intersect_forConcordance.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance)

concordance_bySTR_simpleSTR_v2_allSTR_STR_merge <- read_table("~/Desktop/KU/@research/STR/figure/EH_TRGT_STR_length_compare_byID_merge_allSTR_withconcordance1.txt")
head(concordance_bySTR_simpleSTR_v2_allSTR_STR_merge)
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge %>% select(ID,STR_ID,STR_diff,allele) -> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_foranno

str_anno <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/STR.concordance.INFO.withanno.simpleSTR.txt") %>% 
  mutate(type = ifelse(type == "upstream","promoter",type)) %>%
  mutate(type = ifelse(str_detect(type,"RNA"),"ncRNA",type))
head(str_anno)
str_anno %>% select(STR_ID,type,cyto_type,main_cyto_type,chrom,distance_bycen) -> str_anno
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_foranno %>% left_join(str_anno) -> concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_foranno
concordance_bySTR_simpleSTR_v2_allSTR_STR_merge_foranno %>% group_by(chrom) %>% count(STR_diff != 0) %>% mutate(prop = prop.table(n)) %>%
  #filter(`STR_diff != 0` == "FALSE")
  ggplot(aes(x=chrom,y=n,fill=prop)) + 
  geom_col()

str_anno <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/STR.concordance.INFO.withanno.simpleSTR.txt") %>% 
  mutate(type = ifelse(type == "upstream","promoter",type)) %>%
  mutate(type = ifelse(str_detect(type,"RNA"),"ncRNA",type))

head(str_anno)
str_anno %>% 
  count(type,concordance_range) %>% 
  mutate(type = fct_reorder(type, n, .desc = TRUE)) %>%  # n 기준으로 type을 내림차순 정렬
  ggplot(aes(x=type,y=n,fill=factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill") + 
  theme_step1() + 
  labs(x="Genomic annotation",y=" ") + 
  geom_text(aes(label = scales::comma(n), vjust = ifelse(concordance_range == 1, 0, 0)), 
            position = position_fill(vjust = 0), size = 5) +  # 조건에 labs(x="Genomic annotation",y="Proportion of STR") + 
  theme(legend.position = "none")

str_anno %>% filter(chrom == "chrX") %>% #head()
  mutate(chrXtype = case_when(
    end <= 2699520 ~ "PAR1",
    start >= 154931044 & end <= 156030895 ~ "PAR2",
    start >= 2699521 & end <= 154931043 ~ "nonPAR",
    TRUE ~ "other")) %>%
  mutate(par = ifelse(chrXtype == "nonPAR",chrXtype,"PAR")) %>% #count(par)
  count(type,concordance_range) %>%
  ggplot(aes(x=type,y=n,fill=as.factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill") + 
  theme_step1() + 
  labs(x="Genomic annotation",y=" ") + 
  geom_text(aes(label = scales::comma(n), vjust = ifelse(concordance_range == 1, 0, 0)), 
            position = position_fill(vjust = 0), size = 5) +  # 조건에 labs(x="Genomic annotation",y="Proportion of STR") + 
  theme(legend.position = "none")

str_anno %>% filter(chrom == "chrX") %>% #head()
  mutate(chrXtype = case_when(
    end <= 2699520 ~ "PAR1",
    start >= 154931044 & end <= 156030895 ~ "PAR2",
    start >= 2699521 & end <= 154931043 ~ "nonPAR",
    TRUE ~ "other")) %>%
  mutate(par = ifelse(chrXtype == "nonPAR",chrXtype,"PAR")) %>% #count(par)
  count(par,concordance_range) %>%
  ggplot(aes(x=par,y=n,fill=as.factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill") + 
  theme_step1() + 
  labs(x=" ",y=" ") + 
  geom_text(aes(label = scales::comma(n), vjust = ifelse(concordance_range == 1, 0, 0)), 
            position = position_fill(vjust = 0), size = 5) +  # 조건에 labs(x="Genomic annotation",y="Proportion of STR") + 
  theme(legend.position = "none") 

head(str_anno)
str_anno %>% count(cyto_type,concordance_range) %>% #count(cyto_type)
  ggplot(aes(x=factor(cyto_type,levels=c("gneg","gpos25","gpos50","gpos75","gpos100","gvar","acen","stalk")),y=n,fill=as.factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill") + 
  theme_step1() + 
  labs(x="Cytoband",y=" ") + 
  geom_text(aes(label = scales::comma(n), vjust = ifelse(concordance_range == 1, 0, 0)), 
            position = position_fill(vjust = 0), size = 5) +  # 조건에 labs(x="Genomic annotation",y="Proportion of STR") + 
  theme(legend.position = "none")
str_anno %>% count(concordance_range)
str_anno %>% mutate(concordance_rate_range = case_when(
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
  ggplot(aes(x=concordance_rate_range,y=distance_bycen)) +
  geom_violin()

head(str_anno)
lm(concordance_rate~type+main_cyto_type+cyto_type+distance_bycen+GC+chrom,data = str_anno) -> a
summary(a)$coefficients %>% as.data.frame() -> asso_result
summary(a)$coefficients[2,]
asso_result$type <- row.names(asso_result)
head(asso_result)
colnames(asso_result) <- c("Estimate","Std.Error","t.value","P","type")
asso_result$type
asso_result$`Pr(>|t|)`
#write.table(asso_result,"~/Desktop/KU/@research/STR/figure/figure3/Asso.annovar.annotation.txt",col.names = T,row.names = F,quote = F,sep = "\t")
asso_result %>% filter(type != "(Intercept)") %>% #filter(-log10(`Pr(>|t|)`) > 10 )
  mutate(main_type = case_when(
  str_detect(type,"chrom") ~ "Chromosome",
  str_detect(type,"main_cyto_type") ~ "Cytoband",
  str_detect(type,"cyto_type") ~ "cytoband2",
  str_detect(type,"type") ~ "STR annotation",
  TRUE ~ "Distance to centromere")) %>% 
  mutate(new_type = case_when(
    str_detect(type,"main_cyto_type") ~ str_replace_all(type,"main_cyto_type",""),
    str_detect(type,"cyto_type") ~ str_replace_all(type,"cyto_type",""),
    str_detect(type,"type") ~ str_replace_all(type,"type",""),
    str_detect(type,"chrom") ~ str_replace_all(type,"chrom",""),
    TRUE ~ type)) %>% filter(main_type != "cytoband2",type!="GC") %>% #count(type)
  arrange(main_type) %>% #head()
  ggplot(aes(y = -log10(P), x = factor(main_type,levels=c("Chromosome","Cytoband","STR annotation","Distance to centromere")), color = main_type,size=Estimate)) +  # main_type에 따른 색상
  geom_point() + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  # y=0.5 수평선 추가
  geom_text(aes(label = ifelse(-log10(P) > 4, new_type, "")), size = 4, color = "black", vjust = 1) +  # -log10이 10 이상인 경우에만 텍스트 표시
  #scale_color_manual(values = c("chrom" = "blue", "cytoband" = "green", "Annotation" = "red")) +  # main_type에 따른 색상 지정
  #ylim(c(1,30)) +
  labs(y="-log10(P)") + 
  coord_flip() + 
  theme_step1() + 
  guides(color = "none")  + 
  theme(axis.title.y = element_blank())


ggarrange(p2.1,nrow = 1,labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> f2.up
ggarrange(p2.4,nrow = 1,labels = c('B'), font.label = list(size = 28), label.y = 1.01) -> f2.mid

ggarrange(p2.1,nrow = 1,labels = c('A'), font.label = list(size = 28), label.y = 1.01) %>%
  annotate_figure(top = text_grob("Whole Simple STR", size = 28, face = "bold")) -> f2.up
ggarrange(p2.2,p2.3,nrow = 1,labels = c('C','D'), font.label = list(size = 28), label.y = 1.01,
          widths = c(3.5, 1)) %>% 
  annotate_figure(top = text_grob("Simple STR in Chromosome X", size = 28, face = "bold")) -> f2.down
