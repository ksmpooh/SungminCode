## figure 1
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggbreak)
library(png)
library(grid)
###
#version = 'v1'

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
trgt_pass_prop <- read.table("figure1.trgt_pass_prop.txt",header = T)
eh_pass_prop <- read.table("figure1.eh_pass_prop.txt",header = T)

trgt_pass_prop %>% rbind(eh_pass_prop) %>% #head()
  ggplot(aes(x=type,y=prop,fill=type)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  theme_step1()

#trgt_pass_prop <- read.table("figure1.trgt_pass_prop.txt",header = T)
trgt_pass_prop <- read.table("figure1.trgt_ap.mean_prop.txt",header = T)
eh_pass_prop <- read.table("figure1.eh_pass_prop.txt",header = T)
head(trgt_pass_prop)


trgt_pass_prop %>% 
  ggplot(aes(y=AP_mean,x=type,fill="#F8766D")) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  theme_step1() + 
  #  ylim(c(0.99,1)) + 
  labs(y="mean of AP") + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) -> p2.1

eh_pass_prop %>% #head()
  ggplot(aes(x=type,y=prop,fill=type)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black") + 
  
  #  theme_bw() + 
  theme_step1() +
  # ylim(c(0.99,1)) + 
  labs(y="Proportion of PASS (%)") + 
  scale_fill_manual(values = c("Short-read" = "#619CFF")) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) -> p2.2


ggarrange(p2.1,p2.2) -> p2
p2
#### venn
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_STR_ID_forVenn.txt",header = T)
head(venn_rawdata)
library(VennDiagram) 
venn_data <- list(
  STR_DB = venn_rawdata %>% row.names(),
  trgt = venn_rawdata %>% select(TRGT) %>% na.omit() %>% row.names(),
  eh = venn_rawdata %>% select(EH) %>% na.omit() %>% row.names()
)
head(venn_data)
#venn_data %>% select(TRGT) %>% na.omit() %>% row.names() %>% list()
p3 <- venn.diagram(x = venn_data, filename = NULL,
                            category.names = c("STR DB", "Long-read", "Short-read"),
                            col = c("#440154AA", "#21908DAA", "#FDE725AA"), # 원의 색상을 반투명하게 설정
                            fill = c("#440154AA", "#21908DAA", "#FDE725AA"), 투명하게 설정
                            alpha = 0.3, # 원의 반투명도를 설정
                            cex = 1.9, # 카테고리 레이블 크기
                            cat.cex = 1.9, # 집합 이름의 글자 크기
                            cat.fontface = "bold", # 글자 폰트 스타일 설정
                            cat.fontfamily = "Arial")
cowplot::plot_grid(p3)


#ex_fig1 = ggarrange(ex_fig1_a, ex_fig1_b,ex_fig1_c, nrow = 1, labels = c('a','b','c'), font.label = list(size = 28), label.y = 1.01,
 #                   widths = c(1.2,1.4,0.5))

ggarrange(p2.1,p2.2,p3,nrow = 1, labels = c('b','c','d'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1,1,2)) -> f1_up


###

length_common_motif <- read.table("figure1.length_common_motif.count.txt",header = T)
ru.scale = c(seq(2,14),"15+")
ru.scale
#table(length_common_motif$length)

length_common_motif %>% 
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"15+")) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  labs(x="Length of Repeat Units",y="# of STRs") -> p4.1
p4.1

length_common_motif %>% #filter(is.na(length))
  mutate(RU.length = ifelse(RU.length %in% ru.scale,RU.length,"15+")) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  theme_step1() + 
  scale_y_continuous(labels = scales::comma) + 
  xlim(c(seq(7,14),"15+")) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())-> p4.2
p4.2

p4 <- ggdraw() +
  draw_plot(p4.1, 0, 0, 1, 1) +           # p1을 원본 크기로 배치
  draw_plot(p4.2, 0.4, 0.4, 0.6, 0.6)     # p2를 오른쪽 위에 작게 배치

#combined_plot

ggarrange(p2.1,p2.2,p3,nrow = 1, labels = c('b','c','d'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1,1,2)) -> f1_up
#ggarrange(f1_up,p4,nrow = 2,ncol = , labels = c('b','c','d'), font.label = list(size = 28), label.y = 1.01)
cowplot::plot_grid(f1_up,p4,ncol = 1)



#####
common_motif_count_frequency <- read.table("figure1.frequency_common_motif.count.txt",header = T)





common_motif_count_frequency %>%
  ggplot(aes(x=factor(new_n,c("1","2","3","4","5","6","7","8","9","10","11+")),y=n)) + 
  geom_bar(stat='identity') + 
  labs(x="Frequency of Repeat Units",y="# of STRs") +
  theme_step1() + 
  scale_y_break(c(1000, 4800),ticklabels = c(1,2,3,4900,5000)) + 
    theme(#axis.title.y = element_blank(),
          legend.position = "none") ->  p5.1

p5.1

common_motif_count_frequency %>%
  ggplot(aes(x=factor(new_n,c("1","2","3","4","5","6","7","8","9","10","11+")),y=n)) + 
  geom_bar(stat='identity') + 
  theme_step1() + 
  xlim(c("4","5","6","7","8","9","10")) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) ->  p5.2


p5 <- ggdraw() +
  draw_plot(p5.1, 0, 0, 1, 1) +           # p1을 원본 크기로 배치
  draw_plot(p5.2, 0.4, 0.4, 0.6, 0.6)

p5

ggarrange(p4,p5,nrow = 1, labels = c('e','f'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1.3,1)) -> f1_down

cowplot::plot_grid(f1_up,f1_down,ncol = 1) -> p1_right

############
img <- readPNG("f1.workflow.png")
g <- rasterGrob(img, interpolate=TRUE)
p <- ggplot() +
  theme_void()
p + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  theme_step1() -> p1

ggarrange(p1,nrow = 1, labels = c('a'), font.label = list(size = 28), label.y = 1.01) -> p1
          
###

cowplot::plot_grid(p1,p1_right,ncol = 2,rel_widths = c(3,5)) 

