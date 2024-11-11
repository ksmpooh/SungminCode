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
#version = 'v3'

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
trgt_pass_prop %>% rbind(eh_pass_prop) %>% group_by(type) %>%
  summarise(mean(prop))

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
venn_rawdata %>% count(str_detect(STR_DB,"chr"))
#321253+ 73
#321326
venn_rawdata %>% na.omit() %>% dim() #311542
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
                   category.names = c("STR DB", "LRS", "SRS"),
                   cat.cex = 2,
                   col = c("grey", "#F8766D", "#619CFF"),
                   fill = c("grey", "#F8766D", "#619CFF"), # 원 내부 색상도 반투명하게 설정
                   alpha = 0.3, # 원의 반투명도를 설정
                   cex = 1.5, # 카테고리 레이블 크기름의 글자 크기
                   #cat.fontface = "bold", # 글자 폰트 스타일 설정
                   cat.fontfamily = "Arial")
cowplot::plot_grid(p3)


#ex_fig1 = ggarrange(ex_fig1_a, ex_fig1_b,ex_fig1_c, nrow = 1, labels = c('a','b','c'), font.label = list(size = 28), label.y = 1.01,
#                   widths = c(1.2,1.4,0.5))

ggarrange(p2,p3,nrow = 1, labels = c('B','C'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1,1.5)) -> f1.up
f1.up

###


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

ggarrange(p4,labels = c('D'), font.label = list(size = 28), label.y = 1.01) -> f1.mid
head(length_unique_motif.count.LvsS)
length_unique_motif.count.LvsS %>% group_by(type) %>% filter(RU.length  > 14) %>% summarise(sum(n))
#length_unique_motif.count.LvsS %>% count(RU.length)
## p5 
notinterserct <- read.table("~/Desktop/KU/@research/STR/figure/figure1.notintersect_TRGT_EH_merge_processing.txt",header = T)
head(notinterserct)

notinterserct %>% group_by(g) %>% count(STR_length > 100)

notinterserct %>% mutate(g = ifelse(g == "TRGT","LRS","SRS")) %>%
  ggplot(aes(x=g,y=STR_length,fill=g)) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  theme_step1() + 
  labs(y="Mean of STR length (bp)") +
  theme(legend.position = "none",
        axis.title.x = element_blank()) -> p5
  

notinterserct %>% mutate(GC = GC * 100) %>% 
  mutate(g = ifelse(g == "TRGT","LRS","SRS")) %>%
  ggplot(aes(x=g,y=GC,fill=g)) + 
  geom_violin() + geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.5, color = "black") + 
  scale_fill_manual(values = c("LRS" = "#F8766D","SRS" = "#619CFF")) + 
  theme_step1() + 
  labs(y="GC content (%)") +
  theme(legend.position = "none",
        axis.title.x = element_blank()) -> p6


## add GC hexa
gc_STRcallrateforSRS <- read.table("~/Desktop/KU/@research/STR/figure/figure1.trgt_notintersect_withEHcallrate.txt",header = T)
head(gc_STRcallrateforSRS)


gc_STRcallrateforSRS %>%
  ggplot(aes(x=EH_pass,y=GC)) +
  geom_hex(bins = 15) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "STR call rate for SRS", y = "GC contents (%)", fill = "# of STRs") +  # 축 및 범례 라벨
  theme_step1() -> p6


model <- lm(GC ~ EH_pass,data=gc_STRcallrateforSRS)

# 회귀 계수 확인
summary(model)

# 회귀식 확인
coef(model)







notinterserct %>% mutate(GC = GC * 100) %>% #head()
  filter(GC > 50) %>% count(g)

notinterserct %>% mutate(GC = GC * 100) %>% #head()
  filter(GC > 75) %>% count(g)
#1   EH 159
#2 TRGT 327
notinterserct %>% mutate(GC = GC * 100) %>% #head()
   filter(GC > 50) %>% count(g)
#1   EH 456
#2 TRGT 868


p4
ggarrange(p4,labels = c('D'), font.label = list(size = 28), label.y = 1.01) -> f1.mid
#f1.mid


ggarrange(p5,p6,nrow = 1,labels = c('E','F'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1,1)) -> f1.down

f1.down



#combined_plot

f1.down

cowplot::plot_grid(f1.up,f1.mid,f1.down,ncol = 1,rel_heights = c(1,1,1)) -> f1.right
f1.right



#####

############
pdf_file <- "~/Desktop/KU/@research/STR/figure/mainworkflow.v3.pro.pdf"

pdf_convert(pdf_file, format = "png", pages = 1, dpi = 300, filenames = "mainworkflow.v3.pro_covert.png")
g <- png::readPNG("~/Desktop/KU/@research/STR/figure/mainworkflow.v3.pro_covert.png")
g <- rasterGrob(g, interpolate = TRUE)

#g
p <- ggplot() +
  theme_void()
p + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  theme_step1() -> p1
p1
ggarrange(p1,nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> p1
p1
###
cowplot::plot_grid(p1,f1.right,ncol = 2,rel_widths = c(2.5,5)) 

png("~/Desktop/KU/@research/STR/figure/final/f1.png", width = 2000, height = 1000)
cowplot::plot_grid(p1,f1.right,ncol = 2,rel_widths = c(3,3.5)) 
dev.off()

cowplot::plot_grid(p1,f1.right,ncol = 2,rel_widths = c(3,3.5)) 

ggsave("~/Desktop/KU/@research/STR/figure/final/f1.png",
       dpi=300, dev='png', width=26,height=13, units="in")


 