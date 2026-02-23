## figure 3 APscore bamQC
## simple STR
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)


###
#version = 'v3'

setwd('~/Desktop/KU/@research/STR/figure/figure3/')
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
# A tibble: 1 × 2
  `mean(concordance_rate)` `sd(concordance_rate)`
                     <dbl>                  <dbl>
1                    0.954                0.00378
"

APmean_byID_simpleSTR <- read_table("f3.AP_mean_byID_simpleSTR.afterchrXQC.txt")
head(APmean_byID_simpleSTR)
colnames(APmean_byID_simpleSTR)[2] <- "APmean"

APmean_byID_simpleSTR %>% summarise(mean(APmean),sd(APmean))
"
# A tibble: 1 × 2
  `mean(APmean)` `sd(APmean)`
           <dbl>        <dbl>
1          0.990    0.0000728
"
###
STR.allele.count.match.byAP <- read_table("f3.STR.allele.count.match.byAP.txt")
STR.allele.count.byAP <- read_table("f3.STR.allele.count.byAP.txt")
sample.concordance.rate.byAP <- read_table("f3.sample.concordance.rate.byAP.txt")
sample.concordance.rate.bySTR_length_AP <- read_table("f3.sample.concordance.rate.bySTR_length_AP.txt")
match.count.bySTR_length_AP <- read_table("f3.match.count.bySTR_length_AP.txt")

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
p3.1
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
  mutate(check = ifelse(check == "1","LRS=SRS  ","LRS≠SRS  ")) %>% #head()
  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(TRGT_AP,levels=rev(c("[0~0.25)","[0.25~0.5)","[0.5~0.75)","[0.75~1]"))),y=n,fill=check)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="proportion of LRS vs SRS",x="\nRange of AP score") + 
  scale_fill_manual(values = c("LRS=SRS  " = "#00BFC4","LRS≠SRS  " = "#F8766D")) +  # fill 색 반대로 설정
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
ggarrange(p1.1,p1.2,p1.3,nrow = 1,widths = c(1,0.5,0.4),align = "h") -> p1
p1


##
concordacne_byID_AP0.5_simpleSTR_afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure3/f3.concordacne_byID_AP0.5_simpleSTR_afterchrXQC.txt")
concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC  <- read_table("~/Desktop/KU/@research/STR/figure/figure3/f3.concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC.txt")

head(concordacne_byID_AP0.5_simpleSTR_afterchrXQC)
head(concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC)

head(concordacne_byID_AP0.5_simpleSTR_afterchrXQC)

head(concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC)
concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC %>%
  ggplot(aes(x=factor(AP0.5,levels=c("[0~0.5]","(0.5~1]")),y=concordance_rate,fill=AP0.5)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black",alpha=0.8) +
  facet_grid(~factor(STR_length,levels=c("[0~50)","[50~Inf)"))) +
  ylim(c(0,1)) + 
  labs(y="Concordance",x="Range of AP score",title = "STR Length") +
  theme_step1() + 
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 18,family = 'Arial'),
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5)) -> p2
p2

ggarrange(p1,p2,nrow = 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01, widths = c(1,0.5),align = "h") -> f4.1
f4.1


##########
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC")
head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)



concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>%  
  pivot_longer(trgt_meanmapq:eh_meanmapq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name,alpha=0.8)) + 
  geom_density() + 
  labs(y="Density",x="Mean of MapQ by STR") +
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme_step1() + 
  theme(legend.position = 'none') -> p3.1
p3.1

concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% 
  select(STR_ID,concordance_rate,trgt_meanbaseq,eh_meanbaseq) %>% 
  filter(concordance_rate == 0) %>%
  pivot_longer(trgt_meanbaseq:eh_meanbaseq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name,alpha=0.8)) + 
  geom_density() + 
  labs(y="Density",x="Mean of BaseQ by STR") +
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme_step1() + 
  theme(legend.position = 'none') -> p3.2

concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% 
  select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>% #head()
  pivot_longer(trgt_meanmapq:eh_meanmapq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS  ","SRS  ")) %>%
  ggplot(aes(x=value,fill=name)) + 
  geom_density(alpha=1) + 
  scale_fill_manual(values = c("LRS  " = "#F8766D","SRS  " = "#619CFF")) + 
  theme_step1() + 
  theme(legend.title = element_blank()) -> p3.3

p3.3 <- cowplot::get_legend(p3.3)

ggarrange(p3.1,p3.2,p3.3, nrow= 1, labels = c('C','D',""), font.label = list(size = 28), label.y = 1.01,widths = c(1,1,0.2)) -> f4.2


cowplot::plot_grid(f4.1,f4.2,ncol = 1,rel_heights = c(1,1)) 



png("~/Desktop/KU/@research/STR/figure/final/f3.png", width = 1000, height = 500)
cowplot::plot_grid(f4.1,f4.2,ncol = 1,rel_heights = c(1,1)) 
dev.off()





head(APmean_byID_bySTRtype)
APmean_byID_bySTRtype %>% 
  mutate(new_g = paste0(patho,"\n(",type,")")) %>% select(ID,new_g,APmean) %>%
  rbind(APmean_byID_allSTR %>% mutate(new_g = "Overall") %>% select(ID,new_g,APmean)) %>% #count(new_g)
  ggplot(aes(x=factor(new_g,levels=c("Overall","Normal\n(simpleSTR)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")),y=APmean,fill=factor(new_g,levels=c("Overall","Normal\n(simpleSTR)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")))) +
  geom_violin() +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1.5, color = "black",alpha=0.6) +
  theme_step1() + 
  labs(y="Mean of AP score by Sample",title = "") + 
  theme(legend.position = 'none',
        axis.title.x = element_blank()) -> p1
p1  





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
p2
ggarrange(p1,p2, nrow= 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,widths = c(1,1)) -> f4.1
f4.1


###
head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC")
head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)



concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>%  count()
  pivot_longer(trgt_meanmapq:eh_meanmapq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name,alpha=0.8)) + 
  geom_density() + 
  labs(y="Density",x="Mean of MapQ by STR") +
  theme_step1() + 
  theme(legend.position = 'none')

concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% #select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  #  filter(concordance_rate == 0) %>%
  ggplot(aes(x=trgt_meanbaseq,y=eh_meanbaseq)) + 
  geom_point()


concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>%
  pivot_longer(trgt_meanmapq:eh_meanmapq) %>% mutate(name=ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  ggplot(aes(x=value,fill=name,alpha=0.8)) + 
  geom_density() + 
  labs(y="Density",x="Mean of MapQ by STR") +
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
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
  theme(legend.position = 'none') #-> p4

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

####


STR.allele.count.match.byAP <- read_table("f3.STR.allele.count.match.byAP.txt")
STR.allele.count.byAP <- read_table("f3.STR.allele.count.byAP.txt")
sample.concordance.rate.byAP <- read_table("f3.sample.concordance.rate.byAP.txt")
sample.concordance.rate.bySTR_length_AP <- read_table("f3.sample.concordance.rate.bySTR_length_AP.txt")
match.count.bySTR_length_AP <- read_table("f3.match.count.bySTR_length_AP.txt")

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
# A tibble: 4 × 3
  TRGT_AP     mean      sd
  <chr>      <dbl>   <dbl>
1 [0.25~0.5) 0.355 0.0195 
2 [0.5~0.75) 0.340 0.00933
3 [0.75~1]   0.957 0.00379
4 [0~0.25)   0.177 0.00969
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


