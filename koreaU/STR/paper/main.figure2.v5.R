## figure 2
# + ciros annotation + with comporae eh count vs TRGT
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggbreak)
library(png)
library(grid)
library(ggpie)
library(ggforce)
library(patchwork)
library(viridis)

###
#version = 'v5'

extrafont::font_import(pattern = "Arial", prompt = F)
extrafont::loadfonts()
setwd('~/Desktop/KU/@research/STR/figure/figure2_withchrX/')
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


final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")



#c("#F8766D", '#00BA38', '#619CFF')

### figure 2 
# f1 sample ID
concordance_byID <- read.table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.concordance_byID_simpleSTR.afterchrXQC.txt",header = T)
head(concordance_byID)
concordance_byID %>% summarise(mean(concordance_rate),sd(concordance_rate))
#mean(concordance_rate) sd(concordance_rate)
#1              0.9538486          0.003777398

concordance_byID %>% 
  ggplot(aes(x = "", y = concordance_rate)) + 
  geom_violin(trim=FALSE,fill="gray") + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  labs(y="Concordance",x="") +
  theme_step1() + 
  theme(axis.ticks = element_blank(),
        axis.title.x = element_blank()) -> p1

#rm(eh_trgt_merge_simple_pass_intersect_forConcordance)

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1")
cor(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC$TRGT_STR,eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC$EH_STR)
#[1] 0.6098551
eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX <- read_table("eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX.txt")
eh_trgt_merge_simple_pass_intersect_meanSTRlengthINFO_afterQCchrX <- read_table("eh_trgt_merge_simple_pass_intersect_meanSTRlengthINFO_afterQCchrX.txt")
head(eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX)
dim(eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX)
eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX %>% filter(STR_ID == "RFC1")
head(eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX)


simple_concordance <- read.table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.concordance_bySTR_simpleSTR.afterchrXQC.txt",header = T)
head(simple_concordance)
mean(simple_concordance$concordance_rate) #0.9539438
sd(simple_concordance$concordance_rate) #0.1340218
simple_concordance %>% mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>% group_by(g) %>% #head()
  summarise(mean(concordance_rate),sd(concordance_rate))
"
  g      `mean(concordance_rate)` `sd(concordance_rate)`
  <chr>                     <dbl>                  <dbl>
1 normal                    0.954                  0.134
2 patho                     0.939                  0.176
"
#simple_concordance %>% filter(STR_ID == "RFC1")
simple_concordance %>% mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>% 
  mutate(same = ifelse(concordance_rate == 1,"1","0"))%>% count(g,same) %>% group_by(g) %>% mutate(prop=prop.table(n)*100)
head(simple_concordance)
#simple_concordance %>% filter(STR_ID == "RFC1")
#140617 + 170853 = 311470
62
"
  g      same       n  prop
  <chr>  <chr>  <int> <dbl>
1 normal 0     125139  40.2
2 normal 1     186331  59.8
3 patho  0         36  58.1
4 patho  1         26  41.9
"
simple_concordance %>% mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>%
  mutate(same = ifelse(concordance_rate == 1,"1","0"))%>% count(same)  %>% mutate(prop=prop.table(n)*100)
'
  same      n     prop
1    0 125175 40.18046
2    1 186357 59.81954
'

simple_concordance %>% #mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>%
  mutate(same = ifelse(concordance_rate == 0,"O","X"))%>% count(same)  %>% mutate(prop=prop.table(n)*100)
#concordacen 0
#409 0.1312867

#png("~/Desktop/KU/@research/STR/figure/sup.figure/STR_length.png", width = 1000, height = 200)
#png("~/Desktop/KU/@research/STR/figure/sup.figure/STR_length.png",res = 300)

###### sup
eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX %>% 
  #  mutate(facet = ifelse(mean_TRGT_STR < 150,1,2)) %>% 
  #filter(facet == 2) %>%
  ggplot(aes(x = mean_TRGT_STR, y = mean_EH_STR)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "blue",alpha=0.8, se = FALSE) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  labs(y ="SRS",x="LRS",title = "Mean of STR Repeat Count") +
  theme_step1() + 
  coord_fixed(ratio = 1,ylim = c(0,400)) + 
  theme(legend.position = "none",
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5)) -> p2
#p2

ggarrange(p1,p2, nrow = 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,widths = c(1,5)) -> f2.1
f2.1

head(eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX)
head(eh_trgt_merge_simple_pass_intersect_meanSTRlengthINFO_afterQCchrX)
eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX %>% 
  summarise(max(max_TRGT_STR)) # 3658
eh_trgt_merge_simple_pass_intersect_meanSTRcountINFO_afterQCchrX %>% 
  filter(max_TRGT_STR == 3658) #chr21_10272501_10272513
final_ref %>% filter(ID == "chr21_10272501_10272513") #TA
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(TRGT_STR == 3658)
#1 NIH23F1144740 chr21_10272501_10272513     3658     23    0.56   0.529      2
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "chr21_10272501_10272513",ID == "NIH23F1144740")
#ID            STR_ID                  TRGT_STR EH_STR TRGT_AM TRGT_AP allele
#<chr>         <chr>                      <dbl>  <dbl>   <dbl>   <dbl>  <dbl>
# 1 NIH23F1144740 chr21_10272501_10272513        6      4    0.5    1          1
#2 NIH23F1144740 chr21_10272501_10272513     3658     23    0.56   0.529      2

#chr21	10272501	10272513	16	13	100	16	38.8	60

eh_trgt_merge_simple_pass_intersect_meanSTRlengthINFO_afterQCchrX %>% 
  summarise(max(max_TRGT_STR_legnth)) # 12535
eh_trgt_merge_simple_pass_intersect_meanSTRlengthINFO_afterQCchrX %>% 
  filter(max_TRGT_STR_legnth == 12535) %>% select(STR_ID) #chr22_16317305_16317315
#head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "chr22_16317305_16317315") %>% 
  summarise(max(TRGT_STR))
#2507
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "chr22_16317305_16317315") %>% filter(TRGT_STR == "2507")
#1 NIH23F1602908 chr22_16317305_16317315     2507      2    0.64   0.664      2
final_ref %>% filter(ID == "chr22_16317305_16317315") #ATTCC

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "chr22_16317305_16317315",ID == "NIH23F1602908")
#ID            STR_ID                  TRGT_STR EH_STR TRGT_AM TRGT_AP allele
#<chr>         <chr>                      <dbl>  <dbl>   <dbl>   <dbl>  <dbl>
#  1 NIH23F1602908 chr22_16317305_16317315        2      2   -1      1          1
#2 NIH23F1602908 chr22_16317305_16317315     2507      2    0.64   0.664      2
##rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
#chr22	16317305	16317315	18	11	100	18	39.3	59.7


eh_trgt_merge_simple_pass_intersect_meanSTRlengthINFO_afterQCchrX %>% left_join(simple_concordance) %>%
  filter(concordance_rate != 1) %>%
  ggplot(aes(x = mean_TRGT_STR_legnth, y = mean_EH_STR_legnth)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  #geom_errorbar(aes(ymin = min_EH_STR_legnth, ymax = max_EH_STR_legnth), width = 0,alpha = 0.5) +  # y ???????????? ??????????????????
  #geom_errorbarh(aes(xmin = min_TRGT_STR_legnth, xmax = max_TRGT_STR_legnth), height = 0,alpha = 0.5) +  # x ???????????? ??????????????????
  labs(y ="SRS",x="LRS",title = "Mean of STR Length") +
  theme_step1() + 
  theme(legend.position = "none",
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5))+
  coord_fixed(ratio = 1) -> sp.2.2
#sup_str_length
#cowplot::plot_grid(f2.1,f2.2,f2.3,f2.4,ncol = 1,rel_heights = c(1,1,1,1.8)) 
#dev.off()
#sp.2.2



############



concordanceRange_bySTR_simpleSTR <- read.table("f2.concordance.range.bySTR_simpleSTR.afterchrXQC.txt",header = T)
head(concordanceRange_bySTR_simpleSTR)
head(concordanceRange_bySTR_simpleSTR)
concordanceRange_bySTR_simpleSTR %>%
  mutate(facet = ifelse(concordance_rate_range %in% c("[0.9,1)","1"),2,1)) %>% #head()
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                       "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                       "[0.8,0.9)", "[0.9,1)","1")),y=n, fill = n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  #  scale_fill_viridis_c(option = "viridis") +  
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Concordance") + 
  geom_text(aes(label = scales::comma(n), y = n), 
            size = 5, family = "Arial", vjust = 0) +
  facet_row(~ facet, scales = "free", space = "free") + 
  theme(strip.text = element_blank(),
        #plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'),
        legend.position = "none") -> p3
ggarrange(p3,nrow = 1, labels = c('C'), font.label = list(size = 28), label.y = 1.01) -> f2.2
f2.2

STR.allele.count.match.bySTR_length <- read_table("f2.STR.allele.count.match.bySTR_length.txt")
STR.allele.count.bySTR_length <- read_table("f2.STR.allele.count.bySTR_length.txt")
sample.concordance.rate.bySTR_length <- read_table("f2.sample.concordance.rate.bySTR_length.txt")
sample.concordance.rate.bySTR_length_GC <- read_table("f2.sample.concordance.rate.bySTR_length_GC.txt")
match.count.bySTR_length_GC <- read_table("f2.match.count.bySTR_length_GC.txt")

head(STR.allele.count.match.bySTR_length)
head(STR.allele.count.bySTR_length)
head(sample.concordance.rate.bySTR_length)
head(sample.concordance.rate.bySTR_length_GC)
head(match.count.bySTR_length_GC)


head()
#1.	#F8766D (Red)
#2.	#7CAE00 (Green)
#3.	#00BFC4 (Cyan)
#4.	#C77CFF (Purple)

head(STR.allele.count.match.bySTR_length)
STR.allele.count.match.bySTR_length
'
  STR_length check        n   prop
  <chr>      <dbl>    <dbl>  <dbl>
1 [0~50)         0  1632775 0.0419
2 [0~50)         1 37354282 0.958 
3 [100~150)      0    16505 0.532 
4 [100~150)      1    14508 0.468 
5 [150~Inf)      0    41350 0.980 
6 [150~Inf)      1      863 0.0204'
STR.allele.count.match.bySTR_length %>% group_by(STR_length) %>% #filter(STR_length != "Overall") %>% 
  mutate(check = ifelse(check == "1","LRS=SRS  ","LRS≠SRS  ")) %>% #head()
  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),y=n,fill=check)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="Proportion of LRS vs SRS",x="\nRange of STR length (bp)") + 
  scale_fill_manual(values = c("LRS=SRS  " = "#00BFC4","LRS≠SRS  " = "#F8766D")) +  # fill 색 반대로 설정
  geom_text(aes(label = ifelse(str_detect(check,"="),"",paste0(round(prop,1),"%")), y = 0),
            hjust = 0, size = 5)  +
  geom_text(aes(label = ifelse(str_detect(check,"="),paste0(round(prop,1),"%"),""), y = 1),
            hjust = 1, size = 5)  +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #plot.margin = unit(c(0.25, 0.5, 1, 0.5), 'cm'),
        plot.margin = unit(c(5.5, 10, 5.5, 15), "pt"),
        #plot.margin = unit(c(5.5, -10, 5.5, 5.5), "pt"),
        strip.text = element_blank()) -> p4.1
# unit(c(top, right, bottom, left), units)
STR.allele.count.bySTR_length %>% group_by(STR_length) %>% #filter(STR_length != "Overall") %>% 
  summarise(n=sum(n)) %>%  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),y=n)) + 
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
        strip.text = element_blank()) -> p4.2

STR.allele.count.match.bySTR_length %>% group_by(STR_length) %>% #filter(STR_length != "Overall") %>% 
  mutate(check = ifelse(check == "1","LRS=SRS  ","LRS≠SRS  ")) %>% #head()
  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),y=n,fill=check)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="proportion of LRS vs SRS",x="\nRange of STR length (bp)") + 
  scale_fill_manual(values = c("LRS=SRS  " = "#00BFC4","LRS≠SRS  " = "#F8766D")) +  # fill 색 반대로 설정
  geom_text(aes(label = ifelse(str_detect(check,"="),"",paste0(round(prop,1),"%")), y = 0),
            hjust = 0, size = 5)  +
  geom_text(aes(label = ifelse(str_detect(check,"="),paste0(round(prop,1),"%"),""), y = 1),
            hjust = 1, size = 5)  +
  theme(legend.title = element_blank(),
        #legend.position = "none",
        #plot.margin = unit(c(5.5, 10, 5.5, -10), "pt"),
        plot.margin = unit(c(5.5, -10, 5.5, 5.5), "pt"),
        strip.text = element_blank()) -> p4.3

p4.3 <- cowplot::get_legend(p4.3)

#ggarrange(p2.1,p2.2,nrow = 1, labels = c('B'), font.label = list(size = 28), label.y = 1.01, widths = c(0.9,1),align = "h") -> p2
ggarrange(p4.1,p4.2,p4.3,nrow = 1, widths = c(1,0.6,0.3),align = "h") -> p4
p4



###
sample.concordance.rate.bySTR_length



sample.concordance.rate.bySTR_length %>% ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),
                                                    y=concordance_rate,fill=factor(STR_length,levels=c("[0~50)","[50~100)","[100~150)","[150~Inf)")))) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5,alpha=0.8)  + 
  labs(y="Concordance ",x=" ") + 
  ylim(c(0,1)) + 
  coord_flip() + 
  theme_step1() +
  theme(legend.position = "none",
        #        axis.text.x = element_text(color = "white"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        legend.title = element_text()) -> p5
p5
head(sample.concordance.rate.bySTR_length)
sample.concordance.rate.bySTR_length %>% group_by(STR_length)%>%summarise(mean = mean(concordance_rate),sd = sd(concordance_rate))
"
STR_length   mean      sd
<chr>       <dbl>   <dbl>
  1 [0~50)     0.958  0.00343
2 [100~150)  0.468  0.0221 
3 [150~Inf)  0.0204 0.00482
4 [50~100)   0.860  0.0169 
"
ggarrange(p4,p5,nrow = 1, labels = c('D','E'), font.label = list(size = 28), label.y = 1.01, widths = c(2,1)) -> f2.3
f2.3






####f3
head(sample.concordance.rate.bySTR_length_GC)
head(match.count.bySTR_length_GC)

#concordance.range.byID_length_GC_simpleSTR <- read_table("f2.concordance.range.byID_length_GC_simpleSTR.txt")
head(sample.concordance.rate.bySTR_length_GC)
#sample.concordance.rate.bySTR_length_GC %>% group_by(GC) %>% summarise(mean(concordance_rate))


sample.concordance.rate.bySTR_length_GC %>%
  mutate(GC = case_when(
    GC == "[0~0.25)" ~ "0~25",
    GC == "[0.25~0.5)" ~ "25~50",
    GC == "[0.5~0.75)" ~ "50~75",
    TRUE ~ '75~100')) %>% #count(GC)
  ggplot(aes(x=factor(GC,levels=c("0~25","25~50","50~75","75~100")),
             y=concordance_rate,fill=factor(GC,levels=c("0~25","25~50","50~75","75~100")))) + 
  geom_violin() + 
  labs(x="GC content (%)",y="Concordance",title = "STR Length") + 
  theme(legend.position = 'none') + 
  facet_grid(~factor(STR_length,levels=c("[0~50)","[50~100)","[100~150)","[150~Inf)"))) +
  theme_step1() + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18,family = 'Arial'),
        plot.margin = unit(c(10, 5.5, 5.5, 20), "pt"),
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5)) -> p6
p6
ggarrange(p6,nrow = 1, labels = c('F'), font.label = list(size = 28), label.y = 1.01) -> f2.4


cowplot::plot_grid(f2.1,f2.2,f2.3,f2.4,ncol = 1,rel_heights = c(1,1,1,1)) 


### circos
circos_img <- png::readPNG("circos_plot.png")
g_circos <- rasterGrob(circos_img, interpolate = TRUE)

# Circos plot
p7 <- ggplot() +
  theme_void() +
  annotation_custom(g_circos, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

print(p7)
#ggarrange(p5, nrow = 1, labels = c('E'), font.label = list(size = 28), label.y = 1.01) -> p5


####  annotation


anno.concordance <- read_table("f2.annotation.info.concordance.txt")
head(anno.concordance)
tail(anno.concordance)

head(anno.concordance)
anno.concordance %>% filter(concordance != 1) %>% #filter(type== "Cytoband")
  filter(type %in% c("Chromosome", "Arm", "Cytoband", "Annotation")) %>%
  mutate(type = ifelse(type == "Annotation","Genomic\nfunction",type)) %>%
  group_by(type) %>%
  mutate(
    highest = concordance == max(concordance),
    lowest = concordance == min(concordance)
  ) %>% 
  mutate(highest = ifelse(type == "Arm",NA,highest)) %>%
  ggplot(aes(x = type, y = concordance, fill = type, color = type)) +
  geom_point(position = position_jitter(width = 0.05, height = 0)) + 
  geom_text(
    aes(label = ifelse(highest, name, "")),  
    vjust = 2,  
    size = 5
  ) +
  geom_text(
    aes(label = ifelse(lowest, name, "")), 
    vjust = -1, 
    size = 5
  ) +
  ylim(c(0.8,1)) + 
  theme_step1() + 
  coord_flip() + 
  theme(legend.position = 'none',
        axis.title.y = element_blank()) + 
  labs(y = "Concordance") -> p8

p8
ggarrange(p7,p8, nrow = 1, labels = c('E','F'), font.label = list(size = 28), label.y = 1.01,widths = c(1.4,1)) -> f2.5



cowplot::plot_grid(f2.1,f2.2,f2.3,f2.4,f2.5,ncol = 1,rel_heights = c(1,1,1,1,1.2)) 


ggarrange(p1,p2, nrow = 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,widths = c(1,5)) -> sf2.1
ggarrange(sp.2.2, nrow = 1, labels = c('C'), font.label = list(size = 28), label.y = 1.01) -> sf2.2


ggarrange(p3,nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> f2.1
ggarrange(p4,p5,nrow = 1, labels = c('B','C'), font.label = list(size = 28), label.y = 1.01, widths = c(2,1)) -> f2.2
ggarrange(p6,nrow = 1, labels = c('D'), font.label = list(size = 28), label.y = 1.01) -> f2.3
ggarrange(p7,p8, nrow = 1, labels = c('E','F'), font.label = list(size = 28), label.y = 1.01,widths = c(1.4,1)) -> f2.4

#####p4
#######
png("~/Desktop/KU/@research/STR/figure/final/f2.png", width = 1000, height = 1400)
cowplot::plot_grid(f2.1,f2.2,f2.3,f2.4,ncol = 1,rel_heights = c(1,1,1,1.8)) 
dev.off()


cowplot::plot_grid(f2.1,f2.2,f2.3,f2.4,ncol = 1,rel_heights = c(1.1,1,1,1.8)) 
ggsave("~/Desktop/KU/@research/STR/figure/final/f2.png",
       dpi=300, dev='png', width=13,height=17, units="in")




png("~/Desktop/KU/@research/STR/figure/sup.figure/STR_length.png", width = 1000, height = 500)
cowplot::plot_grid(sf2.1,sf2.2,ncol = 1,rel_heights = c(1,1)) 
dev.off()



