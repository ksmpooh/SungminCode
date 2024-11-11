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

###
#version = 'v1'

extrafont::font_import(pattern = "Arial", prompt = F)
extrafont::loadfonts()
setwd('~/Desktop/KU/@research/STR/figure/figure2/')
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

#c("#F8766D", '#00BA38', '#619CFF')

### figure 2

concordance_byID_simpleSTR <- read.table("f2.concordance.range.byID_simpleSTR.txt",header = T)
#concordance_byID_simpleSTR
concordance_byID_simpleSTR %>% ggplot(aes(x = "", y = concordance_rate)) + 
  geom_violin(trim=FALSE,fill="gray") + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  labs(y="Concordance rate") +
  theme_step1() + 
  theme(axis.ticks = element_blank(),
        axis.title.x = element_blank()) -> p1

head(concordance_byID_simpleSTR)
mean(concordance_byID_simpleSTR$concordance_rate)

p1
#head(concordanceRange_bySTR_simpleSTR)
concordanceRange_bySTR_simpleSTR <- read.table("f2.concordance.range.bySTR_simpleSTR.txt",header = T)
concordanceRange_bySTR_simpleSTR %>%
  mutate(facet = ifelse(concordance_rate_range %in% c("[0.9,1)","1"),2,1)) %>% #head()
  ggplot(aes(x=factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                       "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                       "[0.8,0.9)", "[0.9,1)","1")),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Concordance rate") + 
  geom_text(aes(label = scales::comma(n), y = n), # 바의 위에 텍스트를 배치하기 위해 y를 n보다 약간 크게 설정
            size = 5, family = "Arial", vjust = 0) +
  facet_row(~ facet, scales = "free", space = "free") + 
  theme(strip.text = element_blank()) -> p2

p2
head(concordanceRange_bySTR_simpleSTR)
#          widths = c(1,1.5)
ggarrange(p1,p2,nrow = 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01, widths = c(1,4.5)) -> f2.up

f2.up



#####

ru.scale = c(seq(2,19),"20+")
ru.scale1 = c(seq(2,6))
ru.scale2 = c(seq(7,19),"20+")

concordance_bySTR_simpleSTR <- read.table("f2.concordance_bySTR_simpleSTR.txt",header = T)

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
  geom_text(aes(label = ifelse(concordance_rate_range %in% c("x = 0"),n,""),y = 1,vjust = 0),size = 5) + 
  geom_text(aes(label = ifelse(concordance_rate_range %in% c("0 < x < 1"),scales::comma(n),"")),position = position_stack(vjust = 0.1),size = 5) + 
  
  geom_text(aes(label = ifelse(concordance_rate_range %in% c("x = 1"),scales::comma(n),""),y=0.05,vjust = 0),size = 5) + 
  theme_step1() +
  theme(legend.direction = "horizontal",
        legend.position = "bottom") -> p3

head(concordance_bySTR_simpleSTR)
concordance_bySTR_simpleSTR %>% filter(!(str_detect(STR_ID,'chr'))) %>% summarise(mean(concordance_rate)) #0.8717949
concordance_bySTR_simpleSTR %>% filter((str_detect(STR_ID,'chr'))) %>% summarise(mean(concordance_rate)) #0.9239342

concordance_bySTR_simpleSTR %>% filter(!(str_detect(STR_ID,'chr'))) %>% filter(concordance_rate != 1) %>% 
  summarise(mean(concordance_rate)) #0.8205128

p3
#####



eh_trgt_merge_simple_forConcordance_meanRepeatcount_bySTR <- read_table("f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.txt")
head(eh_trgt_merge_simple_forConcordance_meanRepeatcount_bySTR)

eh_trgt_merge_simple_forConcordance_meanRepeatcount_bySTR %>% select(STR_ID:mean_Whole.length) %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length) %>% count(EH < TRGT)

#EH > TRGT: 77315
#EH < TRGT: 41478

eh_trgt_merge_simple_forConcordance_meanRepeatcount_bySTR %>% select(STR_ID:mean_Whole.length) %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length) %>% mutate(diff = TRGT - EH) %>% count(diff <=  -100)

#TRGT: diff > 200 197. >=100 326
#EH: diff > 100 197. >=100 3

eh_trgt_merge_simple_forConcordance_meanRepeatcount_bySTR %>% select(STR_ID:mean_Whole.length) %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length) %>% mutate(diff = TRGT - EH) %>% filter(diff != 0 ) %>%
  filter(str_detect(STR_ID,"chr")) %>% summarise(mean(diff)) # 1.11

eh_trgt_merge_simple_forConcordance_meanRepeatcount_bySTR %>% select(STR_ID:mean_Whole.length) %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length) %>% mutate(diff = TRGT - EH) %>% filter(diff != 0 ) %>%
  filter(!str_detect(STR_ID,"chr")) %>% summarise(mean(diff)) # 3.62


eh_trgt_merge_simple_forConcordance_meanRepeatcount_bySTR %>% select(STR_ID:q3_Whole.length) %>%
  pivot_wider(names_from = name,values_from = mean_Whole.length:q3_Whole.length) %>% #head(10000) %>% #head()
  mutate(facet = ifelse(mean_Whole.length_TRGT > 200,2,1)) %>% #head()
  mutate(q3_Whole.length_TRGT = ifelse(facet == 1 & q3_Whole.length_TRGT > 210,210,q3_Whole.length_TRGT)) %>%
  mutate(q1_Whole.length_TRGT = ifelse(facet == 2 & q1_Whole.length_TRGT < 200,200,q1_Whole.length_TRGT)) %>% #head
  ggplot(aes(x = mean_Whole.length_TRGT, y = mean_Whole.length_EH)) +
  geom_point() +  # 점 그래프
  geom_errorbar(aes(ymin = q1_Whole.length_EH, ymax = q3_Whole.length_EH), width = 0,alpha = 0.2) +  # y 방향 에러바
  geom_errorbarh(aes(xmin = q1_Whole.length_TRGT, xmax = q3_Whole.length_TRGT), height = 0,alpha = 0.2) +  # x 방향 에러바
  labs(x = "Mean of STR length by LRS", y = "Mean of STR length by SRS") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  theme_step1() + 
  facet_row(~ facet, scales = "free", space = "fixed") +
  theme(strip.text = element_blank()) -> p4



ggarrange(p3,p4,nrow = 2, labels = c('C','D'), font.label = list(size = 28), label.y = 1.01) -> f2.mid


f2.up
f2.mid

png("~/Desktop/KU/@research/STR/figure/final/f2.png", width = 1250, height = 1400)
cowplot::plot_grid(f2.up,f2.mid,ncol = 1,rel_heights = c(1,2)) 
dev.off()
