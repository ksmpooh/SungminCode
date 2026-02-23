## figrue 3 adfter chr X QC after RFC1
# pathogenic
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(cowplot)


###
#version = 'v3'

setwd('~/Desktop/KU/@research/STR/figure/figure3/')
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

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro
head(final_ref)
final_ref_pro
head(final_ref) 

final_ref %>% filter(ID == "RFC1") 
final_ref %>% filter(str_detect(MOTIFS,",")) %>% 
  separate_rows(MOTIFS, sep = ",") %>%
  rbind(final_ref %>% filter(ID == "RFC1")) %>% #filter(ID == "RFC1")
  group_by(chrom, start, end, ID) %>%
  mutate(priority = row_number()) %>%
  ungroup() %>%
  mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX"))) %>%
  arrange(chrom, start) %>%
  rename(STR_ID = ID) %>%
  mutate(new_ID = paste0(STR_ID,"_",MOTIFS)) -> final_ref_complex

final_ref_complex %>% filter(STR_ID == "RFC1") 


head(concordance.byID.nonpatho.patho.simpleSTR.complexSTR)
##### p1
concordance.byID.nonpatho.patho.simpleSTR.complexSTR <- read.table("~/Desktop/KU/@research/STR/figure/figure3/concordance.byID.nonpatho.patho.simpleSTR.complexSTR.afterchrXQC.txt",header = T)
concordance.byID.nonpatho.patho.simpleSTR.complexSTR %>% group_by(g,str_type) %>%
  summarise(mean(concordance_rate))

"
#<chr>          <chr>                           <dbl>
  g          str_type     `mean(concordance_rate)`
  <chr>      <chr>                           <dbl>
1 Normal     (simpleSTR)                     0.954
2 Pathogenic (Overall)                       0.941
3 Pathogenic (complexSTR)                    0.959
4 Pathogenic (simpleSTR)                     0.938
"

#concordance.bySTR.nonpatho.patho.simpleSTR.complexSTR <- read.table("concordance.byID.nonpatho.patho.simpleSTR.complexSTR.afterchrXQC.txt",header = T)
#head(concordance.bySTR.nonpatho.patho.simpleSTR.complexSTR)

#concordance.byID.nonpatho.patho.simpleSTR.complexSTR %>% filter(STR_+I)

concordance.byID.nonpatho.patho.simpleSTR.complexSTR %>% filter(str_type %in% c("(simpleSTR)","(Overall)")) %>%
  ggplot(aes(x=factor(g,levels=c("Normal","Pathogenic")),y=concordance_rate,fill=factor(g,levels=c("Normal","Pathogenic")))) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1.5, color = "black",alpha=0.8) + 
  #ylim(c(0.9)) +
  labs(y="Concordance") + 
  theme_step1() + 
  scale_fill_manual(values = c("Normal" = "#1f77b4", "Pathogenic" = "#d62758")) +
  theme(legend.position = 'none',
        axis.title.x=element_blank()) -> p1
p1

########
final_ref_complex %>% filter(STR_ID == "RFC1")
patho_complex_processing <- read_table("complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")
patho_complex_processing %>% filter(STR_ID == "RFC1")
patho_complex_processing %>% mutate(new_ID = ifelse(new_ID == "RFC1","RFC1_AARRG",new_ID)) %>%
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>%
  group_by(STR_ID, new_ID) %>%
  summarise(concordance_rate = mean(check)) %>%  #filter(STR_ID == "RFC1")
  left_join(final_ref_complex) %>%  #filter(STR_ID == "RFC1")
  mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX"))) %>%
  arrange(chrom, start) %>% 
  select(STR_ID,new_ID,concordance_rate,chrom) #%>% writexl::write_xlsx("~/Desktop/KU/@research/STR/table/complex.STR.each.motif.concordance.xlsx")
"
'''
STR_ID   new_ID                           concordance_rate chrom
   <chr>    <chr>                                       <dbl> <fct>
 1 ATXN7    ATXN7_GCA                                   1     chr3 
 2 ATXN7    ATXN7_GCC                                   1     chr3 
 3 CNBP     CNBP_CA                                     1     chr3 
 4 CNBP     CNBP_CAGA                                   1     chr3 
 5 CNBP     CNBP_CAGG                                   1     chr3 
 6 HTT      HTT_CAG                                     0.923 chr4 
 7 HTT      HTT_CCG                                     1     chr4 
 8 RFC1     RFC1_AARRG                                  0.677 chr4 
 9 FXN      FXN_A                                       1     chr9 
10 FXN      FXN_GAA                                     1     chr9 
11 FRA10AC1 FRA10AC1_CCA                                1     chr10
12 FRA10AC1 FRA10AC1_CCG                                1     chr10
13 ATXN8OS  ATXN8OS_CTA                                 1     chr13
14 ATXN8OS  ATXN8OS_CTG                                 1     chr13
15 NOP56    NOP56_CGCCTG                                1     chr20
16 NOP56    NOP56_GGCCTG                                1     chr20
17 PRNP     PRNP_CCTCAGGGCGGTGGTGGCTGGGGGCAG            1     chr20
18 PRNP     PRNP_CCTCATGGTGGTGGCTGGGGGCAG               1     chr20
19 TMEM185A TMEM185A_CGC                                1     chrX 
20 TMEM185A TMEM185A_CGCCGT                             1     chrX 
'''
"
##
pathogenic.all.allele.strlength.compare.txt <- read_table("~/Desktop/KU/@research/STR/figure/figure3/pathogenic.all.allele.strlength.compare.txt")
pathogenic.all.allele.strlength.compare.txt %>% filter(STR_ID == "RFC1")
head(pathogenic.all.allele.strlength.compare.txt)
pathogenic.all.allele.strlength.compare.txt %>% filter(TRGT_STR_length > 1000)
#≠


pathogenic.all.allele.strlength.compare.txt %>% mutate(check = ifelse(check == 0,"LRS≠SRS","LRS=SRS")) %>% #head()
  mutate(STR_type = ifelse(STR_type == "complex","complex STR","simple STR")) %>% #head()
  mutate(STR_length = case_when(
    TRGT_STR_length < 50 ~ "[0~50)",
    TRGT_STR_length < 100 ~ "[50~100)",
    TRGT_STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) -> pathogenic.all.allele.strlength.compare_addSTRlengthRange

pathogenic.all.allele.strlength.compare_addSTRlengthRange %>% 
  count(STR_length) -> a
head(pathogenic.all.allele.strlength.compare_addSTRlengthRange1)
pathogenic.all.allele.strlength.compare_addSTRlengthRange %>% 
  count(STR_length,check) %>% group_by(STR_length) %>%
  mutate(prop = prop.table(n))->b
b
head(a)
head(b)
b %>% mutate(prop=prop.table(n)*100)  %>% #head()
ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),y=n,fill=check)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="Proportion of LRS vs SRS",x="\nRange of STR length (bp)") + 
  scale_fill_manual(values = c("LRS=SRS" = "#00BFC4","LRS≠SRS" = "#F8766D")) +  # fill �� �ݴ�� ����
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
p2.1
head(a)
a %>% group_by(STR_length) %>% #filter(STR_length != "Overall") %>% 
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
        strip.text = element_blank()) -> p2.2


b %>%# head()
  mutate(check = ifelse(check == "LRS=SRS","LRS=SRS  ","LRS≠SRS  ")) %>% #head()
  ggplot(aes(x=factor(STR_length,levels=rev(c("[0~50)","[50~100)","[100~150)","[150~Inf)"))),y=n,fill=check)) + 
  geom_bar(stat = "identity",position = "fill") + 
  coord_flip() + 
  theme_step1() + 
  labs(y="Proportion of LRS vs SRS",x="\nRange of STR length (bp)") + 
  scale_fill_manual(values = c("LRS=SRS  " = "#00BFC4","LRS≠SRS  " = "#F8766D")) +  # fill �� �ݴ�� ����
  geom_text(aes(label = ifelse(str_detect(check,"="),"",paste0(round(prop,2),"%")), y = 0),
            hjust = 0, size = 5)  +
  geom_text(aes(label = ifelse(str_detect(check,"="),paste0(round(prop,2),"%"),""), y = 1),
            hjust = 1, size = 5)  +
  theme(legend.title = element_blank(),
        #legend.position = "none",
        #plot.margin = unit(c(0.25, 0.5, 1, 0.5), 'cm'),
        plot.margin = unit(c(5.5, 10, 5.5, 15), "pt"),
        #plot.margin = unit(c(5.5, -10, 5.5, 5.5), "pt"),
        strip.text = element_blank()) -> p2.3

p2.3 <- cowplot::get_legend(p2.3)

#ggarrange(p2.1,p2.2,nrow = 1, labels = c('B'), font.label = list(size = 28), label.y = 1.01, widths = c(0.9,1),align = "h") -> p2
ggarrange(p2.1,p2.2,p2.3,nrow = 1, widths = c(1,0.6,0.3),align = "h") -> p2
p2





pathogenic.all.allele.strlength.compare.txt %>% 
  group_by(STR_type,STR_ID) %>%
  summarise(concordance = mean(check)) %>%
  ggplot(aes(x=STR_type,y=concordance,fill=STR_type)) + 
  geom_violin()








pathogenic.all.allele.strlength.compare.txt %>% mutate(check = ifelse(check == 0,"LRS≠SRS","LRS=SRS")) %>%
  mutate(STR_type = ifelse(STR_type == "complex","complex STR","simple STR")) %>%
  ggplot(aes(x=TRGT_STR_length,y=EH_STR_length,shape=STR_type,color = check)) +
  geom_point(alpha=0.7,size=5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  labs(x="LRS",y="SRS",title = "STR Length by Allele") + 
  scale_color_manual(values = c("LRS≠SRS" = "#a30000", "LRS=SRS" = "#00BFC4")) +  # 색상 설정
  theme_step1() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(family = 'Arial', size = 16, color = 'black',hjust = 0.5),
        legend.direction = "horizontal", 
        legend.spacing.x = unit(1, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        legend.box = "horizontal") + 
  coord_equal(ratio = 1) -> p2

p2

patho_complex_processing %>% filter(STR_ID == "HTT")

patho_complex_processing %>% 
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>%
  group_by(STR_ID, new_ID) %>%
  summarise(concordance_rate = mean(check)) %>% 
  left_join(final_ref_complex %>%   mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX")))) %>% # head()
  group_by(STR_ID) %>% 
  mutate(MOTIFS = factor(MOTIFS, levels = unique(final_ref_complex$MOTIFS[order(final_ref_complex$priority)]))) %>% 
  filter(STR_ID != "PRNP") %>%
  ggplot(aes(x = MOTIFS, y = concordance_rate, fill = STR_ID)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~factor(STR_ID,levels = c("ATXN7","CNBP","HTT","FXN","FRA10AC1","ATXN8OS","NOP56",#"PRNP",
                                       "TMEM185A")),scales = "free_x",ncol = 4) + 
  theme_step1() + 
  labs(y="Concordance") + 
  theme(legend.position = 'none',
        axis.title.x = element_blank()) -> p3


ggarrange(p1,p2, nrow= 2, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,heights = c(1.4,2))




#### mmapQ


#png("~/Desktop/KU/@research/STR/figure/final/f3.png", width = 800, height = 600)
png("~/Desktop/KU/@research/STR/figure/final/f3.png", width = 8, height = 6, units = "in", res = 300)
ggarrange(p1,p2, nrow= 2, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,heights = c(1.4,2))
dev.off()


ggarrange(p1,p2, nrow= 2, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,heights = c(1,1.2))# -> f3

ggsave("~/Desktop/KU/@research/STR/figure/final/f3.png",
       dpi=300, dev='png', height=8, width=8, units="in")

ggsave("~/Desktop/KU/@research/STR/figure/final/f3.png", f3,bg="transparent",
       dpi=300, dev='png', height=8, width=8)

ggsave("~/Desktop/KU/@research/STR/figure/final/f3.pdf", f3,
       dpi=300, dev='pdf', height=8, width=8)



ggarrange(p1,p2, nrow= 1, labels = c('A','B'), font.label = list(size = 28), label.y = 1.01,widths = c(1,1.2)) -> f3_up
f3_up

ggarrange(f3_up,f3_mid,ncol = 1, font.label = list(size = 28), label.y = 1.01,heights = c(1.1,1.5))
