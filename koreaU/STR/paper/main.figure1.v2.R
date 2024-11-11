## figure 1
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggbreak)
library(png)
library(grid)
library(ggpie)
library(ggforce)
###
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

#c("#F8766D", '#00BA38', '#619CFF')

### figure 1

setwd("~/Desktop/KU/@research/STR/figure/")
trgt_pass_prop <- read.table("~/Desktop/KU/@research/STR/figure/figure1.trgt_pass_prop.v2.txt",header = T)
##write.table(eh_ID_accu,"~/Desktop/KU/@research/STR/figure/figure1.eh_pass_prop.v2.txt",col.names = T,row.names = F,quote = F,sep = "\t")
eh_pass_prop <- read.table("~/Desktop/KU/@research/STR/figure/figure1.eh_pass_prop_onlySpan.v2.txt",header = T)

trgt_pass_prop$type = "LRS"
eh_pass_prop$type = "SRS"


#head(trgt_pass_prop)
trgt_pass_prop %>% rbind(eh_pass_prop) %>% #head()
  ggplot(aes(x=type,y=prop,fill=type)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  labs(y="Proportion of PASS (%)") + 
  theme_step1() + 
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank()) -> p2
         
p2




#### venn
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
head(venn_rawdata)
head(venn_rawdata)

venn_rawdata %>% dim() #321326
venn_rawdata %>% select(TRGT) %>% na.omit() %>% dim() # 315839
venn_rawdata %>% select(EH) %>% na.omit() %>% dim() #314244
 
library(VennDiagram) 
venn_data <- list(
  STR_DB = venn_rawdata %>% row.names(),
  trgt = venn_rawdata %>% select(TRGT) %>% na.omit() %>% row.names(),
  eh = venn_rawdata %>% select(EH) %>% na.omit() %>% row.names()
)
head(venn_data)
#venn_data %>% select(TRGT) %>% na.omit() %>% row.names() %>% list()
p3 <- venn.diagram(x = venn_data, filename = NULL,
                   category.names = c("STR DB", "LRS", "SRS"),
                   col = c("grey", "#F8766D", "#619CFF"), # 원의 색상을 반투명하게 설정
                   fill = c("grey", "#F8766D", "#619CFF"), # 원 내부 색상도 반투명하게 설정
                   alpha = 0.3, # 원의 반투명도를 설정
                   cex = 1.9, # 카테고리 레이블 크기
                   cat.cex = 1.9, # 집합 이름의 글자 크기
                   cat.fontface = "bold", # 글자 폰트 스타일 설정
                   cat.fontfamily = "Arial")
cowplot::plot_grid(p3)


#ex_fig1 = ggarrange(ex_fig1_a, ex_fig1_b,ex_fig1_c, nrow = 1, labels = c('a','b','c'), font.label = list(size = 28), label.y = 1.01,
#                   widths = c(1.2,1.4,0.5))

ggarrange(p2,p3,nrow = 1, labels = c('B','C'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1,1.5)) -> f1_up
f1_up

###

length_common_motif <- read.table("figure1.length_common_motif.count.txt",header = T)
ru.scale = c(seq(2,19),"20+")
ru.scale1 = c(seq(2,6))
ru.scale2 = c(seq(7,19),"20+")

sum(length_common_motif$n)
head(length_common_motif)
length_common_motif %>% mutate(prop = prop.table(n)*100) %>% filter(!(RU.length %in% c(1:14))) %>%
  summarise(sum(n))
#length_common_motif %>% is.na()
length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(facet = ifelse(RU.length %in% ru.scale1,1,2)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Length of Repeat Units") + 
  facet_row(~ facet, scales = "free", space = "free") + 
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme(strip.text = element_blank())

#dev.off()
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
  theme(strip.text = element_blank()) -> p4



p4

#plot_grid(f1_up,f1_mid,ncol = 1)

#########
ru.scale = c(seq(2,19),"20+")
ru.scale1 = c(seq(2,6))
ru.scale2 = c(seq(7,19),"20+")



library(ggforce)

length_unique_motif.count.LvsS<- read.table("~/Desktop/KU/@research/STR/figure/figure1.length_unique_motif.count.LvsS.txt",header = T)
length_unique_motif.count.LvsS %>% #head()
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"20+")) %>% 
  mutate(facet = ifelse(RU.length %in% ru.scale1,1,2)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n,fill=type)) + 
  geom_bar(stat = 'identity',position = 'dodge') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(y="# of STRs",x="Length of Repeat Units") + 
  facet_row(~ facet, scales = "free", space = "free") + 
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme(legend.position = "none",
    strip.text = element_blank()) -> p4


## 

p4
p5  
### freqeucny repeat unit p5
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
  theme(legend.position = 'none',
    strip.text = element_blank()) -> p5

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
  theme(strip.text = element_blank())


ggarrange(p4,p5,ncol = 1,labels = c('D','E'), font.label = list(size = 28), label.y = 1.01,
          heights = c(1, 1)) -> f1.mid

f1.mid


  
#combined_plot




######################## pie

annot_common_str <-read.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.CommonSTR.forpie.txt",header = T)
annot_onlyTRGT_str<- read.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.onlyTRGT_STR.forpie.txt",header = T)
annot_onlyEH_str<- read.table("~/Desktop/KU/@research/STR/figure/figrue1.annotation.onlyEH_STR.forpie.txt",header = T)


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



cowplot::plot_grid(f1_up,f1.mid,f1.down,ncol = 1,rel_heights = c(1,2,1)) -> f1_right
f1_right



#####

############
img <- readPNG("~/Desktop/KU/@research/STR/figure/f1.workflow.png")
g <- rasterGrob(img, interpolate=TRUE)
p <- ggplot() +
  theme_void()
p + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  theme_step1() -> p1

ggarrange(NULL,nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> p1

###
f1_right
cowplot::plot_grid(p1,f1_right,ncol = 2,rel_widths = c(3,5)) 

