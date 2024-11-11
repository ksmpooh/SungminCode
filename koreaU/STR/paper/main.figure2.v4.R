## figure 2
# + ciros annotation
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
#version = 'v3'

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

#c("#F8766D", '#00BA38', '#619CFF')

### figure 2
simple_concordance <- read.table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.concordance_bySTR_simpleSTR.afterchrXQC.txt",header = T)
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



concordanceRange_bySTR_simpleSTR <- read.table("f2.concordance.range.bySTR_simpleSTR.afterchrXQC.txt",header = T)
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
        legend.position = "none") -> p1
ggarrange(p1,nrow = 1, labels = c('A'), font.label = list(size = 28), label.y = 1.01) -> f2.1
f2.1

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
      strip.text = element_blank()) -> p2.1
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
        strip.text = element_blank()) -> p2.2

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
        strip.text = element_blank()) -> p2.3

p2.3 <- cowplot::get_legend(p2.3)

#ggarrange(p2.1,p2.2,nrow = 1, labels = c('B'), font.label = list(size = 28), label.y = 1.01, widths = c(0.9,1),align = "h") -> p2
ggarrange(p2.1,p2.2,p2.3,nrow = 1, widths = c(1,0.6,0.3),align = "h") -> p2
p2


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
        legend.title = element_text()) -> p3
p3
head(sample.concordance.rate.bySTR_length)
sample.concordance.rate.bySTR_length %>% group_by(STR_length)%>%summarise(mean = mean(concordance_rate),sd = sd(concordance_rate))
#STR_length   mean      sd
#<chr>       <dbl>   <dbl>
  #1 [0~50)     0.929  0.00342
#2 [100~150)  0.445  0.0224 
#3 [150~Inf)  0.0163 0.00441
#4 [50~100)   0.801  0.0168 



ggarrange(p2,p3,nrow = 1, labels = c('B','C'), font.label = list(size = 28), label.y = 1.01, widths = c(2,1)) -> f2.2
f2.2


ggarrange(f2.up,f2.mid,ncol = 1, font.label = list(size = 28), label.y = 1.01, heights  = c(1,1))


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
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5)) -> p4
p4
ggarrange(p4,nrow = 1, labels = c('D'), font.label = list(size = 28), label.y = 1.01) -> f2.3


cowplot::plot_grid(f2.1,f2.2,f2.3,ncol = 1,rel_heights = c(1,1,1)) 


### circos
circos_img <- png::readPNG("circos_plot.png")
g_circos <- rasterGrob(circos_img, interpolate = TRUE)

# Circos plot
p5 <- ggplot() +
  theme_void() +
  annotation_custom(g_circos, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

print(p5)
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
  labs(y = "Concordance") -> p6

p6
ggarrange(p5,p6, nrow = 1, labels = c('E','F'), font.label = list(size = 28), label.y = 1.01,widths = c(1.4,1)) -> f2.4



cowplot::plot_grid(f2.1,f2.2,f2.3,f2.4,ncol = 1,rel_heights = c(1,1,1,1.2)) 


#####p4
#######
png("~/Desktop/KU/@research/STR/figure/final/f2.png", width = 1000, height = 1300)
cowplot::plot_grid(f2.1,f2.2,f2.3,f2.4,ncol = 1,rel_heights = c(1,1,1,1.8)) 
dev.off()

