### sple or extend data
library(tidyverse)
library(data.table)
library(ggbreak)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(grid)
library(ggpie)
library(patchwork)

#version = 'v2'

extrafont::font_import(pattern = "Arial", prompt = F)
extrafont::loadfonts()

#setwd('~/Dropbox/SMC_AD_WGS_paper/Data/SMC_cwas_results_20240416/')
# Set colors

# Function definition
theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        legend.background = element_rect(color = 'black'), legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        axis.title.y = element_text(margin = margin(r = 10, unit = "pt"))) }




## eh type distribution
setwd("~/Desktop/KU/@research/STR/figure/sup.figure/")


### mapping dataeph
head(depth)
depth <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/mapping.depth.txt") 
depth %>% group_by(platfrom) %>%
  summarise(mean(check),sd(check))
'
platfrom `mean(check)` `sd(check)`
<chr>            <dbl>       <dbl>
1 LRS               33.0        3.68
2 SRS               40.9        2.89
'
read_table("~/Desktop/KU/@research/STR/figure/sup.figure/mapping.depth.txt") %>%
  ggplot(aes(x=platfrom,y=check,fill=platfrom))+ 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black",alpha=0.8) + 
  labs(y="Mean of mapping depth (X)") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x = element_blank())


## eh type distribution

eh_qc_type_count <- read.table("~/Desktop/KU/@research/STR/figure/sup.figure/eh_qc_alleletype_count.bysample.txt",header = T)



eh_qc_type_count %>% mutate(type = ifelse(type == "./.","Missing value",type)) %>%
  ggplot(aes(x=ID,y=n,fill=type)) + 
  geom_bar(position="stack", stat="identity")

table(df_filter_withtype$type)

eh_qc_type_count %>% mutate(type = ifelse(type == "./.","Missing value",type)) %>% filter(FILTER == "LowDepth") %>%
  ggplot(aes(x=ID,y=n,fill=type)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of LowDepth STR") + 
  guides(fill=guide_legend(title="Read type")) + 
  #  scale_fill_manual(values = c("Missing value" = "grey")) + 
  theme(axis.text.x = element_blank())

library(unikn)
#pal_unikn
palette(usecol(pal_unikn, n = 9))    

palette(rainbow(n = 9)) -> a
c("grey",rev(a)) -> a
a

unique(eh_qc_type_count$type)[2:10] %>% sort() -> b
c("Missing",b) -> b
b
a
eh_qc_type_count %>% mutate(type = ifelse(type == "./.","Missing",type)) %>% filter(FILTER == "LowDepth") %>%
  ggplot(aes(x=ID,y=n,fill=factor(type, levels = b))) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of STRs (LowDepth)") + 
  scale_fill_manual(values = a) + 
  theme_step1() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 20)),
        legend.position = 'none') -> p1.2

p1.2



eh_qc_type_count %>% mutate(type = ifelse(type == "./.","Missing",type)) %>% filter(FILTER != "LowDepth") %>%
  ggplot(aes(x=ID,y=n,fill=factor(type, levels = b))) + 
  geom_bar(position="stack", stat="identity") + 
  guides(fill=guide_legend(title="Read type")) + 
  scale_fill_manual(values = a[2:10]) + 
  theme_step1() + 
  scale_y_break(c(30000, 317000),scales = c(1,1)) + 
  xlab("Sample") +
  ylab("# of STRs (PASS)") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5,margin = margin(r = 10)),
        legend.position = 'none') -> p1.1
p1.1  

eh_qc_type_count %>% mutate(type = ifelse(type == "./.","Missing",type)) %>% filter(FILTER == "LowDepth") %>%
  ggplot(aes(x=ID,y=n,fill=factor(type, levels = b))) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of STRs (LowDepth)") + 
  scale_fill_manual(values = a) +
  theme_step1() + 
  guides(fill=guide_legend(title="Read Type by EH")) + 
  theme(
    legend.title = element_text(hjust = 0.5)) -> p1.3

p1.legend <- get_legend(p1.3)


#ggarrange(p1.1,p1.2,ncol = 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01) -> p1.left
p1.left <- p1.1 / p1.2 + plot_layout(heights = c(1, 1))
ggarrange(p1.left,nrow = 1,
          #widths = c(3.5,1),
          common.legend = TRUE, 
          legend.grob = p1.legend,
          legend = "right") -> p1
p1
ggsave("p1.plot.jpg", plot = p1, width = 10, height = 7, dpi = 600)



### annotation
annot_common_str <-read.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.CommonSTR.forpie.txt",header = T)
annot_onlyTRGT_str<- read.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.onlyTRGT_STR.forpie.txt",header = T)
annot_onlyEH_str<- read.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.onlyEH_STR.forpie.txt",header = T)

head(annot_onlyTRGT_str)
ggdonut(annot_common_str,group_key = 'type',count_type = 'full',
        label_info = c("group","ratio"), label_type = "horizon",
        label_size = 4, label_pos = "in", label_threshold = 10) + 
  ggtitle("Common STR") + # 제목 추가
  theme(legend.position = 'none',plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(hjust = 0.5, vjust = 1,family = 'Arial', size = 18, color = 'black',margin = margin(2, 0, 2, 0))) -> p6.1



# 최종 결과 출력

ggdonut(annot_onlyTRGT_str,group_key = 'type',count_type = 'full',
        label_info = c("group","ratio"), label_type = "horizon",
        label_size = 4, label_pos = "in", label_threshold = 10) + 
  ggtitle("Only in LRS") + # 제목 추가
  theme(legend.position = 'none',plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(hjust = 0.5, vjust = 1,family = 'Arial', size = 18, color = 'black',margin = margin(2, 0, 2, 0))) -> p6.2


ggdonut(annot_onlyEH_str,group_key = 'type',count_type = 'full',
        label_info = c("group","ratio"), label_type = "horizon",
        label_size = 4, label_pos = "in", label_threshold = 10) + 
  ggtitle("Only in SRS") + # 제목 추가
  theme(legend.position = 'none',plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(hjust = 0.5, vjust = 1,family = 'Arial', size = 18, color = 'black',margin = margin(2, 0, 2, 0))) -> p6.3

ggarrange(p6.1,p6.2,p6.3,nrow = 1,labels = c('F'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1, 1, 1)) -> f1.down

ggarrange(p6.1,p6.2,p6.3,nrow = 1,labels = c('A','B','C'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1, 1, 1)) -> f1.down

f1.down

#### concordance by mean length



#### concordacen by annotation INFO, Simple STR
str_anno <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/STR.concordance.INFO.withanno.simpleSTR.txt") %>% 
  mutate(type = ifelse(type == "upstream","promoter",type)) %>%
  mutate(type = ifelse(str_detect(type,"RNA"),"ncRNA",type))


str_anno %>% 
  count(type,concordance_range) %>% 
  mutate(type = fct_reorder(type, n, .desc = TRUE)) %>%  # n 기준으로 type을 내림차순 정렬
  ggplot(aes(x=type,y=n,fill=factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill") + 
  theme_step1() + 
  labs(x="Genomic annotation",y=" ") + 
  geom_text(aes(label = scales::comma(n), vjust = ifelse(concordance_range == 1, 0, 0)), 
            position = position_fill(vjust = 0), size = 5) +  # 조건에 labs(x="Genomic annotation",y="Proportion of STR") + 
  theme(legend.position = "none") -> p2.1
  
str_anno %>% filter(chrom == "chrX") %>% #head()
   count(concordance_rate == 1) %>%mutate(prop = prop.table(n))
#1 FALSE                    5318 0.661
#2 TRUE                     2729 0.339

str_anno %>% filter(chrom != "chrX") %>% #head()
  count(concordance_rate == 1) %>%mutate(prop = prop.table(n))
#<lgl>                    <int> <dbl>
#  1 FALSE                   136013 0.448
#2 TRUE                    167474 0.552

str_anno %>% filter(chrom == "chrX") %>% #head()
  mutate(chrXtype = case_when(
    end <= 2699520 ~ "PAR",
    start >= 154931044 & end <= 156030895 ~ "PAR",
    start >= 2699521 & end <= 154931043 ~ "nonPAR",
    TRUE ~ "other")) %>% group_by(chrXtype) %>% count(concordance_rate == 1) %>%mutate(prop = prop.table(n))
#
#chrXtype `concordance_rate == 1`     n  prop
#<chr>    <lgl>                   <int> <dbl>
#1 PAR      FALSE                      54 0.443
#2 PAR      TRUE                       68 0.557
#3 nonPAR   FALSE                    5264 0.664
#4 nonPAR   TRUE                     2661 0.336

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
  theme(legend.position = "none") -> p2.3

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
  theme(legend.position = "none") -> p2.2

str_anno %>% count(cyto_type,concordance_range) %>% #count(cyto_type)
  ggplot(aes(x=factor(cyto_type,levels=c("gneg","gpos25","gpos50","gpos75","gpos100","gvar","acen","stalk")),y=n,fill=as.factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill") + 
  theme_step1() + 
  labs(x="Cytoband",y=" ") + 
  geom_text(aes(label = scales::comma(n), vjust = ifelse(concordance_range == 1, 0, 0)), 
            position = position_fill(vjust = 0), size = 5) +  # 조건에 labs(x="Genomic annotation",y="Proportion of STR") + 
  theme(legend.position = "none") -> p2.4

ggarrange(p2.1,nrow = 1,labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> f2.up
ggarrange(p2.4,nrow = 1,labels = c('B'), font.label = list(size = 28), label.y = 1.01) -> f2.mid

ggarrange(p2.1,nrow = 1,labels = c('A'), font.label = list(size = 28), label.y = 1.01) %>%
  annotate_figure(top = text_grob("Whole Simple STR", size = 28, face = "bold")) -> f2.up
ggarrange(p2.2,p2.3,nrow = 1,labels = c('C','D'), font.label = list(size = 28), label.y = 1.01,
          widths = c(3.5, 1)) %>% 
  annotate_figure(top = text_grob("Simple STR in Chromosome X", size = 28, face = "bold")) -> f2.down



cowplot::plot_grid(f2.up,f2.mid,f2.down,ncol = 1,rel_heights = c(1,1,1)) 





# good 
b_chrX %>% mutate(par = ifelse(chrXtype == "nonPAR",chrXtype,"PAR")) %>% #count(par)
  count(par,concordance_range) %>%
  ggplot(aes(x=par,y=n,fill=as.factor(concordance_range))) +
  geom_bar(stat = "identity",position = "fill")


# 고려
b_chrX %>% mutate(par = ifelse(chrXtype == "nonPAR",chrXtype,"PAR")) %>% #count(par)
  count(par,type,concordance_range) %>% #head()
  ggplot(aes(x=type,y=n,fill=concordance_range)) +
  geom_bar(stat = "identity",position = "fill") + 
  facet_grid(~par,scales = "free")

head(b)

b %>% mutate(concordance_range = ifelse(concordance_rate != 1,0,1)) %>% 
  count(type,concordance_range) %>% 
  ggplot(aes(x=type,y=n,fill=concordance_range)) +
  geom_bar(stat = "identity",position = "fill")

#### bam Q
concordance_bySTR_simpleSTR.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC")
head(concordance_bySTR_simpleSTR.afterchrXQC)
concordance_bySTR_simpleSTR.afterchrXQC %>%
  pivot_longer(cols = c("trgt_meanmapq","eh_meanmapq")) %>% mutate(name = ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name)) + 
  geom_density(alpha =0.5) + 
  facet_wrap(~factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                      "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                      "[0.8,0.9)", "[0.9,1)","1")),scales = 'free_y',nrow = 4) + 
  labs(x="Mean of Mapping Quality by STR") + 
  theme_bw()

head(concordance_bySTR_simpleSTR.afterchrXQC)
concordance_bySTR_simpleSTR.afterchrXQC %>% pivot_longer(cols = trgt_meandepth:eh_meanmapq) %>%  #head()
  mutate(g = ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%  #head()
  group_by(name,g,concordance_rate_range) %>%
  summarise(mean = mean(value)) %>% 
  pivot_wider(names_from = concordance_rate_range,values_from = mean) %>%
  mutate(bamQC = str_split_fixed(name,"_",2)[,2]) %>% rename(platform = g) %>% ungroup() %>%
  arrange(bamQC) %>% select(-name) %>%
  select("bamQC","platform","0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
         "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
         "[0.8,0.9)", "[0.9,1)","1") %>% writexl::write_xlsx("~/Desktop/KU/@research/STR/table/simpleSTR.oribmaQC.info.byconcordacne.xlsx")
  #summarise(mean(trgt_meandepth),mean(trgt_meanbaseq),mean(trgt_meanbaseq),mean(trgt_meanbaseq),mean(trgt_meanbaseq),mean(trgt_meanbaseq),mean(trgt_meanbaseq))
  



concordance_bySTR_simpleSTR.afterchrXQC %>% #head()
  ggplot(aes(x=trgt_meanbaseq,y=eh_meanbaseq)) +
  geom_point(alpha=0.5) + 
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # 회귀선 추가 (선형 회귀)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  ylim(c(0,60)) + xlim(c(0,60)) + 
  labs(x="mean baseQ (LRS)",y="mean baseQ (SRS)") + 
  facet_wrap(~factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                      "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                      "[0.8,0.9)", "[0.9,1)","1")),nrow = 4) + 
  theme_bw()

concordance_bySTR_simpleSTR.afterchrXQC %>% 
  ggplot(aes(x=trgt_maenmapq,y=eh_maenmapq)) +
  geom_point(alpha=0.5) + 
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # 회귀선 추가 (선형 회귀)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  ylim(c(0,60)) + 
  labs(x="mean mapQ (LRS)",y="mean mapQ (SRS)") + 
  facet_wrap(~factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                      "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                      "[0.8,0.9)", "[0.9,1)","1")),nrow = 4)


merge_coverage %>% filter(concordance_rate == 0) %>%
  ggplot(aes(x=EH_maenmapq,y=EH_maenmapq)) +
  geom_point() + 
  labs(x="mean mapQ (LRS)",y="mean mapQ (SRS)") 

    
## complex STR matched
patho_complex_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")

patho_complex_processing %>% #head()
  mutate(check= ifelse(TRGT_STR != EH_STR, 1, 0)) %>% 
  #count(ID,STR_ID)
  group_by(ID,STR_ID,allele) %>%
  summarise(check_rate = sum(check)) %>% #filter(STR_ID == "CNBP")
  group_by(STR_ID) %>% filter(STR_ID != "PRNP") %>%
  count(check_rate) %>% #head()
  ggplot(aes(x=factor(STR_ID,levels = c("ATXN7","HTT","FXN","FRA10AC1","ATXN8OS","NOP56",#"PRNP",
                                        "TMEM185A","CNBP")),y=n,fill=rev(factor(check_rate)))) + 
  geom_bar(stat = 'identity',position = "fill")



############# APscore
setwd('~/Desktop/KU/@research/STR/figure/figure4/')
STR.allele.count.match.byAP <- read_table("f4.STR.allele.count.match.byAP.txt")
STR.allele.count.byAP <- read_table("f4.STR.allele.count.byAP.txt")
sample.concordance.rate.byAP <- read_table("f4.sample.concordance.rate.byAP.txt")
sample.concordance.rate.bySTR_length_AP <- read_table("f4.sample.concordance.rate.bySTR_length_AP.txt")
match.count.bySTR_length_AP <- read_table("f4.match.count.bySTR_length_AP.txt")

head(STR.allele.count.match.byAP)
head(STR.allele.count.byAP)
head(sample.concordance.rate.byAP)
head(sample.concordance.rate.bySTR_length_AP)
head(match.count.bySTR_length_AP)


head()
#1.	#F8766D (Red)
#2.	#7CAE00 (Green)
#3.	#00BFC4 (Cyan)
#4.	#C77CFF (Purple)
table(STR.allele.count.match.byAP$TRGT_AP)
head(STR.allele.count.match.byAP)
STR.allele.count.match.byAP %>% group_by(TRGT_AP) %>% #filter(STR_length != "Overall") %>% 
  mutate(check = ifelse(check == "1","LRS=SRS","LRS!=SRS")) %>% #head()
  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(TRGT_AP,levels=rev(c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]"))),y=n,fill=check)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="Proportion of LRS vs SRS",x="\nRange of AP score") + 
  scale_fill_manual(values = c("LRS=SRS" = "#00BFC4","LRS!=SRS" = "#F8766D")) +  # fill ?????? ??????????? ????????????
  geom_text(aes(label = ifelse(str_detect(check,"LRS!=SRS"),"",paste0(round(prop,1),"%")), y = 0),
            hjust = 0, size = 5)  +
  geom_text(aes(label = ifelse(str_detect(check,"LRS!=SRS"),paste0(round(prop,1),"%"),""), y = 1),
            hjust = 1, size = 5)  +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #plot.margin = unit(c(0.25, 0.5, 1, 0.5), 'cm'),
        plot.margin = unit(c(5.5, 10, 5.5, 15), "pt"),
        #plot.margin = unit(c(5.5, -10, 5.5, 5.5), "pt"),
        strip.text = element_blank()) -> p1.1
p1.1
# unit(c(top, right, bottom, left), units)
STR.allele.count.byAP %>% group_by(TRGT_AP) %>% #filter(STR_length != "Overall") %>% 
  summarise(n=sum(n)) %>%  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(TRGT_AP,levels=rev(c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]"))),y=n)) + 
  geom_col(fill = "grey") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="# of STR alleles",x="\n") + 
  geom_text(aes(label = paste0(scales::comma(n),"(",round(prop,1),"%)"), y = 0),
            hjust = 0, size = 5)  +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        #plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        plot.margin = unit(c(5.5, -10, 5.5, -10), "pt"),
        #
        strip.text = element_blank()) -> p1.2

STR.allele.count.match.byAP %>% group_by(TRGT_AP) %>% #filter(STR_length != "Overall") %>% 
  mutate(check = ifelse(check == "1","LRS=SRS","LRS!=SRS")) %>% #head()
  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(TRGT_AP,levels=rev(c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]"))),y=n,fill=check)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="proportion of LRS vs SRS",x="\nRange of AP score") + 
  scale_fill_manual(values = c("LRS=SRS" = "#00BFC4","LRS!=SRS" = "#F8766D")) +  # fill ?????? ??????????? ????????????
  geom_text(aes(label = ifelse(str_detect(check,"!="),"",paste0(round(prop,1),"%")), y = 0),
            hjust = 0, size = 5)  +
  geom_text(aes(label = ifelse(str_detect(check,"!="),paste0(round(prop,1),"%"),""), y = 1),
            hjust = 1, size = 5)  +
  theme(legend.title = element_blank(),
        #legend.position = "none",
        #plot.margin = unit(c(5.5, 10, 5.5, -10), "pt"),
        plot.margin = unit(c(5.5, -10, 5.5, 5.5), "pt"),
        strip.text = element_blank()) -> p1.3

p1.3 <- cowplot::get_legend(p1.3)

#ggarrange(p1.1,p1.2,nrow = 1, labels = c('B'), font.label = list(size = 28), label.y = 1.01, widths = c(1,0.8),align = "h") -> p1
ggarrange(p1.1,p1.2,p1.3,nrow = 1,widths = c(1,0.6,0.3),align = "h") -> p1
p1


###
sample.concordance.rate.bySTR_length



sample.concordance.rate.byAP %>% 
  ggplot(aes(x=factor(TRGT_AP,levels=rev(c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]"))),
             y=concordance_rate,fill=factor(TRGT_AP,levels=rev(c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]"))))) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5,alpha=0.8)  + 
  labs(y="Concordance ",x=" ") + 
  ylim(c(0,1)) + 
  coord_flip() + 
  theme_step1() +
  theme(legend.position = "none",
        #        axis.text.x = element_text(color = "white"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        legend.title = element_text()) -> p2
p2
head(sample.concordance.rate.bySTR_length)
sample.concordance.rate.byAP %>% group_by(TRGT_AP)%>%summarise(mean = mean(concordance_rate),sd = sd(concordance_rate))
"
  TRGT_AP     mean      sd
  <chr>      <dbl>   <dbl>
1 [0.25~0.5) 0.337 0.0182 
2 [0.5~0.75) 0.330 0.00896
3 [0.75~1]   0.927 0.00377
4 [0~0.25)   0.168 0.00980
"


ggarrange(p1,p2,nrow = 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01, widths = c(2,1)) -> f4.1
f4.1

#ggarrange(f2.up,f2.mid,ncol = 1, font.label = list(size = 28), label.y = 1.01, heights  = c(1,1))


####f3
head(sample.concordance.rate.bySTR_length_AP)
head(match.count.bySTR_length_AP)

#concordance.range.byID_length_GC_simpleSTR <- read_table("f2.concordance.range.byID_length_GC_simpleSTR.txt")
head(sample.concordance.rate.bySTR_length_AP)
#sample.concordance.rate.bySTR_length_GC %>% group_by(GC) %>% summarise(mean(concordance_rate))


sample.concordance.rate.bySTR_length_AP %>%
  ggplot(aes(x=factor(TRGT_AP,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")),
             y=concordance_rate,fill=factor(TRGT_AP,levels=c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]")))) + 
  geom_violin() + 
  labs(x="AP score",y="Concordance",title = "STR Length") + 
  theme(legend.position = 'none') + 
  facet_grid(~factor(STR_length,levels=c("[0~50)","[50~100)","[100~150)","[150~Inf)"))) +
  ylim(c(0,1)) + 
  theme_step1() + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18,family = 'Arial'),
        plot.margin = unit(c(10, 5.5, 5.5, 20), "pt"),
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5)) -> p3
p3
ggarrange(p3,nrow = 1, labels = c('C'), font.label = list(size = 28), label.y = 1.01) -> f4.2


cowplot::plot_grid(f4.1,f4.2,ncol = 1,rel_heights = c(1,1)) 

