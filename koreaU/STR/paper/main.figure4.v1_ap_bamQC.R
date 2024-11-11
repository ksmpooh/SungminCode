## figure 4 APscore
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(cowplot)


###
#version = 'v1'

setwd('~/Desktop/KU/@research/STR/figure/figure4/')
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

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro

#1.	#F8766D (Red)
#2.	#7CAE00 (Green)
#3.	#00BFC4 (Cyan)
#4.	#C77CFF (Purple)


# 
#concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC")

## AP
concordance_byID <- read_table("../figure2_withchrX/concordance_rate_byID.afterchrXQC.txt")
concordance_byID
concordance_byID %>% #ungroup() %>%
  summarise(mean(concordance_rate),sd(concordance_rate))
"
  `mean(concordance_rate)` `sd(concordance_rate)`
                     <dbl>                  <dbl>
1                    0.924                0.00376 "

AP_mean_byID_simpleSTR <- read_table("f4.AP_mean_byID_simpleSTR.afterchrXQC.txt")
head(AP_mean_byID_simpleSTR)

AP_mean_byID_simpleSTR %>% summarise(mean(AP_mean))
#`mean(AP_mean)`
#<dbl>
#  1           0.990
AP_mean_byID_simpleSTR %>%
  ggplot(aes(x="",y=AP_mean)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1.5, color = "black",alpha=0.8) +
  labs(y="AP score") +
  theme_step1() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x =  element_blank()) -> p1
  

concordacne_byID_AP0.5_simpleSTR_afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_simpleSTR_afterchrXQC.txt")
concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC  <- read_table("~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC.txt")

head(concordacne_byID_AP0.5_simpleSTR_afterchrXQC)
head(concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC)

head(concordacne_byID_AP0.5_simpleSTR_afterchrXQC)


concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC %>%
  ggplot(aes(x=factor(AP0.5,levels=c("[0~0.5]","(0.5~1]")),y=concordance_rate,fill=AP0.5)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black",alpha=0.8) +
  facet_grid(~factor(STR_length,levels=c("[0~100)","[100~Inf)"))) +
  ylim(c(0,1)) + 
  labs(y="Concordance",x="Range of AP score",title = "STR Length") +
  theme_step1() + 
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 18,family = 'Arial'),
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5)) -> p2

ggarrange(p1,p2, nrow= 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,widths = c(1.5,3)) -> f4.1
f4.1


###
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC")
head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>%
  pivot_longer(trgt_meanmapq:eh_meanmapq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name,alpha=0.8)) + 
  geom_density() + 
  labs(y="Density",x="Mean of MapQ by STR") +
  theme_step1() + 
  theme(legend.position = 'none') -> p3
  

concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% 
  select(STR_ID,concordance_rate,trgt_meanbaseq,eh_meanbaseq) %>% 
  filter(concordance_rate == 0) %>%
  pivot_longer(trgt_meanbaseq:eh_meanbaseq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name,alpha=0.8)) + 
  geom_density() + 
  labs(y="Density",x="Mean of BaseQ by STR") +
  theme_step1() + 
  theme(legend.position = 'none') -> p4

concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% 
  select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>% #head()
  pivot_longer(trgt_meanmapq:eh_meanmapq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name)) + 
  geom_density(alpha=0.8) + 
  theme_step1() + 
  theme(legend.title = element_blank()) -> p5
p5

p5 <- cowplot::get_legend(p5)

ggarrange(p3,p4,p5, nrow= 1, labels = c('C','D',""), font.label = list(size = 28), label.y = 1.01,widths = c(1,1,0.2)) -> f4.2
f4.2

cowplot::plot_grid(f4.1,f4.2,ncol = 1,rel_heights = c(1,1)) 
