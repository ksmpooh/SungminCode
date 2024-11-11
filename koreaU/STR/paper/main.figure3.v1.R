## figrue 3
# circos plot + anno
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
#version = 'v1'

setwd('~/Desktop/KU/@research/STR/figure/figure3/')
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





## p2
#pdf_file <- "~/Desktop/KU/@research/STR/figure/figure2/circos.portrait_pro2.pdf"

#pdf_convert(pdf_file, format = "png", pages = 1, dpi = 300, filenames = "circos.portrait_covert.png")
circos_img <- png::readPNG("~/Desktop/KU/@research/STR/figure/figure2/circos.portrait_covert.png")
g_circos <- rasterGrob(circos_img, interpolate = TRUE)

# Circos plot을 ggplot 객체로 저장
p1 <- ggplot() +
  theme_void() +
  annotation_custom(g_circos, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

print(p1)
ggarrange(p1, nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01)

# 플롯 표시
print(p)





### p3
library(colorspace)

asso_result <- read.table("~/Desktop/KU/@research/STR/figure/figure3/Asso.annovar.annotation.txt",header = T,sep = "\t")
head(asso_result)

asso_result %>% filter(type != "(Intercept)") %>% #filter(-log10(`Pr(>|t|)`) > 10 )
  mutate(main_type = case_when(
    str_detect(type,"chrom") ~ "Chromosome",
    str_detect(type,"main_cyto_type") ~ "Cytoband",
    str_detect(type,"cyto_type") ~ "cytoband2",
    str_detect(type,"type") ~ "STR annotation",
    TRUE ~ "Distance to centromere")) %>% 
  mutate(new_type = case_when(
    str_detect(type,"main_cyto_type") ~ str_replace_all(type,"main_cyto_type",""),
    str_detect(type,"cyto_type") ~ str_replace_all(type,"cyto_type",""),
    str_detect(type,"type") ~ str_replace_all(type,"type",""),
    str_detect(type,"chrom") ~ str_replace_all(type,"chrom",""),
    TRUE ~ type)) %>% filter(main_type != "cytoband2",type!="GC") %>% filter(-log10(P) > 5)
  arrange(main_type) %>% #head()
  ggplot(aes(y = -log10(P), x = factor(main_type,levels=c("Distance to centromere","STR annotation","Chromosome","Cytoband")), color = main_type,size=Estimate)) +  # main_type에 따른 색상
  geom_point() + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  
  geom_text(aes(label = ifelse(-log10(P) > 4, new_type, "")), size = 4, color = "black", vjust = 1) +  
  scale_color_brewer(palette = "Set2") + 
  labs(y="-log10(P)") + 
  coord_flip() + 
  theme_step1() + 
  guides(color = "none")  + 
  theme(axis.title.y = element_blank()) -> p2

ggarrange(p1,p2, nrow= 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,widths = c(1,1.2)) -> f3_up
f3_up

ggarrange(f3_up,f3_mid,ncol = 1, font.label = list(size = 28), label.y = 1.01,heights = c(1.1,1.5))
