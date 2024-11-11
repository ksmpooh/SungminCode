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
#version = 'v2'

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

####f3

concordance_bySTR_wholelengthmean.noconcordance1 <- read.table("~/Desktop/KU/@research/STR/figure/figure2/f2.eh_trgt_merge_simple_pass_intersect_forConcordance_wholelength.bySTR.mean.notconrconrdacne1.v2.txt",header = T)

concordance_bySTR_wholelengthmean.noconcordance1 %>%
ggplot(aes(x= mean_Whole.length_TRGT,y=mean_Whole.length_EH)) + 
  geom_hex(bins = 30) +  # hexbin plot
  scale_fill_gradient(low = "skyblue", high = "orange") +  # 색상 그라데이션
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y 선 추가
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # 회귀선 추가
  labs(x = "Mean of STR length by LRS", y = "Mean of STR length by SRS", fill = "# of STRs") +  # 축 및 범례 라벨
  coord_fixed(ratio = 1) +
  theme_step1() -> p3

#####p4


img <- readPNG("~/Desktop/KU/@research/STR/figure/figure2/circos.test.png")
g <- rasterGrob(img, interpolate=TRUE)
#g
p <- ggplot() +
  theme_void()
p + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  theme_step1() -> p1
p1



ggarrange(p1,nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> p1
p1
###
cowplot::plot_grid(p1,f1.right,ncol = 2,rel_widths = c(3,5)) 



img <- readPNG("~/Desktop/KU/@research/STR/figure/figure2/circos.test.png")
g <- rasterGrob(img, interpolate = TRUE)

# Plot with adjusted boundaries to reduce whitespace
p <- ggplot() +
  theme_void() +
  annotation_custom(g, xmin = -1, xmax = 1, ymin = -1, ymax = 1) +  # Adjust these values as needed
  theme_step1()

# Adjust the plot to remove excess whitespace
p1 <- p + coord_fixed(ratio = 1)  # Maintain aspect ratio, or you can crop
p1
# Combine the plot with ggarrange
ggarrange(p1, nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01)



####
# PDF 파일 경로
library(pdftools)
library(grid)



# PNG 파일 불러오기

### 이거 굿
pdf_file <- "~/Desktop/KU/@research/STR/figure/figure2/circos.portrait_pro2.pdf"

pdf_convert(pdf_file, format = "png", pages = 1, dpi = 300, filenames = "circos.portrait_covert.png")
circos_img <- png::readPNG("~/Desktop/KU/@research/STR/figure/figure2/circos.portrait_covert.png")
g_circos <- rasterGrob(circos_img, interpolate = TRUE)

# Circos plot을 ggplot 객체로 저장
p1 <- ggplot() +
  theme_void() +
  annotation_custom(g_circos, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

print(p1)
ggarrange(p1, nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01,label.x = 0.02)

# 플롯 표시
print(p)
### 이거 굿



f2.up
f2.mid
#######
png("~/Desktop/KU/@research/STR/figure/final/f2.png", width = 1250, height = 1400)
cowplot::plot_grid(f2.up,f2.mid,ncol = 1,rel_heights = c(1,2)) 
dev.off()
