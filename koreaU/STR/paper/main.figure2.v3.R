## figure 2
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
library(moonBook)

###
#version = 'v3'

extrafont::font_import(pattern = "Arial", prompt = F)
extrafont::loadfonts()
setwd('~/Desktop/KU/@research/STR/figure/figure2/')
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

#c("#F8766D", '#00BA38', '#619CFF')

### figure 2
simple_concordance <- read.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt",header = T)
mean(simple_concordance$concordance_rate) #0.9239236
simple_concordance %>% mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>% group_by(g) %>% #head()
  summarise(mean(concordance_rate)) 
# normal                    0.924
# patho                     0.872
simple_concordance %>% mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>% 
  mutate(same = ifelse(concordance_rate == 1,"1","0"))%>% count(g,same) %>% group_by(g) %>% mutate(prop=prop.table(n)*100)
#1 normal 311470
#2  patho     63
#   g      same       n  prop
#1 normal 0     141285  45.4
#2 normal 1     170185  54.6
#3 patho  0         45  71.4
#4 patho  1         18  28.6



concordanceRange_bySTR_simpleSTR <- read.table("f2.concordance.range.bySTR_simpleSTR.txt",header = T)
head(concordanceRange_bySTR_simpleSTR)
concordanceRange_bySTR_simpleSTR %>%
  mutate(facet = ifelse(concordance_rate_range %in% c("[0.9,1)","1"),2,1)) %>% #head()
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                       "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                       "[0.8,0.9)", "[0.9,1)","1")),y=n, fill = n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  #  scale_fill_viridis_c(option = "viridis") +  # viridis 팔레트 적용
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Concordance") + 
  geom_text(aes(label = scales::comma(n), y = n), # 바의 위에 텍스트를 배치하기 위해 y를 n보다 약간 크게 설정
            size = 5, family = "Arial", vjust = 0) +
  facet_row(~ facet, scales = "free", space = "free") + 
  theme(strip.text = element_blank(),
        legend.position = "none") -> p0


concordance_bySTR_wholelengthmean.noconcordance1 <- read.table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.notconrconrdacne1.v2.txt",header = T)
concordance_bySTR_wholelengthmean <- read.table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.txt",header = T)

head(concordance_bySTR_wholelengthmean)
concordance_bySTR_wholelengthmean %>% group_by(name) %>% summarise(max(mean_Whole.length))
concordance_bySTR_wholelengthmean %>% group_by(name) %>% filter(mean_Whole.length == max(mean_Whole.length))
#STR_ID name mean_Whole.length q1_Whole.length q3_Whole.length max_Whole.length min_Whole.length
#1 chr18_74617039_74617087 TRGT          2590.277            2624            2624             4208              224
concordance_bySTR_wholelengthmean %>% filter(STR_ID == "chr18_74617039_74617087")# %>% filter(mean_Whole.length == max(mean_Whole.length))
concordance_bySTR_wholelengthmean %>% group_by(name) %>% filter(max_Whole.length == max(max_Whole.length))
#2 chr22_16317305_16317315 TRGT              384.               10              10            12535               10
# ATTCC
concordance_bySTR_wholelengthmean %>% filter(STR_ID == "chr22_16317305_16317315")# %>% filter(mean_Whole.length == max(mean_Whole.length))


concordance_bySTR_wholelengthmean.noconcordance1 %>%
  ggplot(aes(x= mean_Whole.length_TRGT,y=mean_Whole.length_EH)) + 
  geom_hex(bins = 30) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "Mean of STR length by LRS (bp)", y = "\nMean of STR length\nby SRS (bp)", fill = "# of STRs") +  # 축 및 범례 라벨
  coord_fixed(ratio = 1) +
  theme_step1() -> p1
#ggarrange(p1,nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> f2.up
ggarrange(p0,p1,ncol = 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01) -> f2.up

f2.up
lm(mean_Whole.length_TRGT ~ mean_Whole.length_EH,concordance_bySTR_wholelengthmean.noconcordance1) -> a
summary(a)
head(concordance_bySTR_wholelengthmean.noconcordance1)
concordance_bySTR_wholelengthmean.noconcordance1 %>% count(mean_Whole.length_TRGT >= mean_Whole.length_EH + 100) #326
concordance_bySTR_wholelengthmean.noconcordance1 %>% count(mean_Whole.length_TRGT >= mean_Whole.length_EH + 50) #540

concordance_bySTR_wholelengthmean.noconcordance1 %>%count(mean_Whole.length_TRGT > mean_Whole.length_EH) #41478
concordance_bySTR_wholelengthmean.noconcordance1 %>%count(mean_Whole.length_TRGT < mean_Whole.length_EH) #77315
concordance_bySTR_wholelengthmean.noconcordance1 %>%count(mean_Whole.length_TRGT == mean_Whole.length_EH) #22537
concordance_bySTR_wholelengthmean.noconcordance1 %>%count(mean_Whole.length_TRGT > 100) #575
concordance_bySTR_wholelengthmean.noconcordance1 %>%count(mean_Whole.length_EH > 100) #202
concordance_bySTR_wholelengthmean.noconcordance1 %>%count(mean_Whole.length_TRGT < 100) #575


#####

STR.allele.count.match.bySTR_length %>% filter(STR_length != "Overall") %>% #summarise(n=sum(n))
  group_by(STR_length) %>% summarise(n=sum(n)) %>% mutate(prop=prop.table(n)*100) 
#  STR_length        n    prop
#<chr>         <dbl>   <dbl>
#1 [0~50)     39215348 96.8   
#2 [100~150)     31083  0.0767
#3 [150~Inf)     42271  0.104 
#4 [50~100)    1210588  2.99  
# 40499290= 311533 *130
#head(STR.allele.count.match.bySTR_length)


STR.allele.count.match.bySTR_length %>% group_by(STR_length) %>% filter(STR_length != "Overall") %>% 
  summarise(n=sum(n)) %>%  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),y=n)) + 
  geom_col(fill = "grey") + 
  coord_flip() + 
  #facet_wrap(~factor(STR_length,levels=c("[0~50)","[50~100)","[100~150)","[150~Inf)")), ncol = 1, scales = "free_y") + 
  theme_step1() + 
  labs(y="# of STRs",x="\nRange of STR length (bp)") + 
  geom_text(aes(label = paste0(scales::comma(n),"(",round(prop,1),"%)"), y = 0),
            hjust = 0, size = 5)  +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.margin = unit(c(5.5, -10, 5.5, 5.5), "pt"),
        strip.text = element_blank()) -> p2.1


STR.allele.count.match.bySTR_length %>% group_by(STR_length) %>% filter(STR_length != "Overall") %>% 
  mutate(match = ifelse(match == "TRUE","LRS=SRS  ","LRS≠SRS  ")) %>%
  mutate(prop=prop.table(n)*100)  %>% #head()
  ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),y=n,fill=match)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="Proportion of LRS vs SRS",x="\n") + 
  geom_text(aes(label = ifelse(str_detect(match,"="),"",paste0(round(prop,1),"%")), y = 0),
            hjust = 0, size = 5)  +
  geom_text(aes(label = ifelse(str_detect(match,"="),paste0(round(prop,1),"%"),""), y = 1),
            hjust = 1, size = 5)  +
  theme(legend.title = element_blank(),
    #legend.position = "none",
        #axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
  plot.margin = unit(c(5.5, 10, 5.5, -10), "pt"),
        strip.text = element_blank()) -> p2.2

#ggarrange(p2.1,p2.2,nrow = 1, labels = c('B'), font.label = list(size = 28), label.y = 1.01, widths = c(0.9,1),align = "h") -> p2
ggarrange(p2.1,p2.2,nrow = 1, widths = c(0.9,1),align = "h") -> p2
p2

concordacen_byID <- read_table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance.range.byID_simpleSTR.v2.txt")  %>% filter(STR_length != "Overall")
concordacen_byID %>% ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),
                                y=prop,fill=factor(STR_length,levels=c("[0~50)","[50~100)","[100~150)","[150~Inf)")))) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5)  + 
  labs(y="Concordance ",x=" ") + 
  ylim(c(0,1)) + 
  coord_flip() + 
  theme_step1() +
  theme(legend.position = "none",
        axis.te.x = element_text(color = "white"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        
        legend.title = element_text()) -> p3
p3
head(concordacen_byID)
concordacen_byID %>% group_by(STR_length)%>%summarise(mean = mean(prop),sd = sd(prop))
#<chr>             <dbl>
#1 [0~50)           0.929 
#2 [100~150)        0.445 
#3 [150~Inf)        0.0163
#4 [50~100)         0.801 



ggarrange(p2,p3,nrow = 1, labels = c('C','D'), font.label = list(size = 28), label.y = 1.01, widths = c(2,1)) -> f2.mid
f2.mid


ggarrange(f2.up,f2.mid,ncol = 1, font.label = list(size = 28), label.y = 1.01, heights  = c(1,1))


####f3

concordance.range.byID_length_GC_simpleSTR <- read_table("f2.concordance.range.byID_length_GC_simpleSTR.txt")
head(concordance.range.byID_length_GC_simpleSTR)
concordance.range.byID_length_GC_simpleSTR %>% filter(STR_length != "Overall") %>%
  mutate(GC = case_when(
    GC == "[0~0.25)" ~ "0~25",
    GC == "[0.25~0.5)" ~ "25~50",
    GC == "[0.5~0.75)" ~ "50~75",
    TRUE ~ '75~100')) %>% group_by(GC) %>% summarise(mean(prop))
  

concordance.range.byID_length_GC_simpleSTR %>% filter(STR_length != "Overall") %>%
  mutate(GC = case_when(
    GC == "[0~0.25)" ~ "0~25",
    GC == "[0.25~0.5)" ~ "25~50",
    GC == "[0.5~0.75)" ~ "50~75",
    TRUE ~ '75~100')) %>% #count(GC)
  ggplot(aes(x=factor(GC,levels=c("0~25","25~50","50~75","75~100")),
             y=prop,fill=factor(GC,levels=c("0~25","25~50","50~75","75~100")))) + 
  geom_violin() + 
  labs(x="GC content (%)",y="Concordance",title = "STR Length") + 
  theme(legend.position = 'none') + 
  facet_grid(~factor(STR_length,levels=c("[0~50)","[50~100)","[100~150)","[150~Inf)"))) +
  theme_step1() + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # y=0.5 수평선 추가
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18,family = 'Arial'),
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5)) -> p4

ggarrange(p4,nrow = 1, labels = c('E'), font.label = list(size = 28), label.y = 1.01) -> f2.down
f2.down
cowplot::plot_grid(f2.up,f2.mid,"",f2.down,ncol = 1,rel_heights = c(2,1,1,1)) 




#####p4
#######
png("~/Desktop/KU/@research/STR/figure/final/f2.png", width = 1000, height = 1300,)
cowplot::plot_grid(f2.up,f2.mid,"",f2.down,ncol = 1,rel_heights = c(2,0.9,0.1,1)) 
dev.off()

ggsave("~/Desktop/KU/@research/STR/figure/final/f2.png",
       cowplot::plot_grid(f2.up,f2.mid,f2.down,ncol = 1,rel_heights = c(1,0.8,1)), 
       width = 8, height = 6, dpi = 300)
