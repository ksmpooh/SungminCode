library(tidyverse)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggbreak)
library(png)
library(grid)
###
setwd("~/Desktop/KU/@research/STR/figure/")
df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReaptnumber.with_dfiff.txt")
str_match_score_freqeuncy <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt")
str_RU_count<- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.RUcount.txt",header =T)
str_qc_pass_common <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.txt",header = T)
ref <- read.table("STR_RU.length_GC.info.txt",header = T)
head(str_RU_count)
head(str_qc_pass_common)
head(str_match_score_freqeuncy)
head(df)

head(df)
complex_str <- read.table("~/Desktop/KU/@research/STR/figure/complex_STR.txt",header = T)
dup_str <- read.table("~/Desktop/KU/@research/STR/figure/rm_dup_str_list.txt",header = T)

head(df)
df %>% select(STR) %>% unique() %>% dim()
df %>% filter(!(STR %in% dup_str$STR)) %>% select(STR) %>% unique() %>% dim() #302656
df %>% filter(!(STR %in% dup_str$STR)) %>%
  mutate(match = ifelse((STR1 == 0 & STR2 == 0),"Match",ifelse((STR1 != 0 & STR2 != 0),"Missmatch","One and One"))) %>% #-> df_match
  count(ID,match)-> df_match
head(df_match)
#df_match %>% write.table("figure2.df_match.txt",col.names = T,row.names = F,quote = F,sep = "\t")
df_match %>% mutate(score = ifelse(match == "Match",n*2,ifelse(match == "Missmatch",n*1,0))) %>% #head() 
  group_by(ID) %>% summarise(sum_s = sum(score)) %>% mutate(ID2='ID') %>%
  mutate(score = sum_s/(302656*2)*100) %>% #head()
  ggplot(aes(y=score,x=ID2)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 1.5, color = "black") + 
  labs(y="Proportion of matching STRs (%)") + 
  theme_step1() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

df_match %>% mutate(pct = n/302656) -> a

302656*2

df_match %>% mutate(score = ifelse(match == "Match",n*2,ifelse(match == "Missmatch",n*1,0))) %>% #head() 
  group_by(ID) %>% summarise(sum_s = sum(score)) %>% 
  mutate(score = sum_s/(302656*2)*100) -> a
  

df %>% filter(!(STR %in% dup_str$STR)) %>%
  mutate(match = ifelse((STR1 == 0 & STR2 == 0),"Match",ifelse((STR1 != 0 & STR2 != 0),"Missmatch","One and One"))) %>% #-> df_match
  count(STR,match)-> df_match_str

head(df_match_str)
df_match_str %>% filter(match == "Match")%>% filter(n == 66) %>% dim() #148513 
148513/302656 # 49.1%
df_match_str %>% filter(match == "Missmatch")%>% filter(n == 66) %>% dim() #189
189/302656 # 0.06%

df %>% filter(!(STR %in% dup_str$STR)) %>% mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(EH_STR = (EH_STR1 + EH_STR2)/2*RU.length) %>%
  mutate(TRGT_STR = (TRGT_STR1 + TRGT_STR2)/2*RU.length) %>%
  select(STR,MOTIFS,ID,EH_STR,TRGT_STR) %>% group_by(STR) %>% 
  summarise(EH = mean(EH_STR),TRGT = mean(TRGT_STR)) -> df_str_length

head(df_str_length)
#df_str_length %>% write.table("figure2.df_str_length_longvsshort.txt",col.names = T,row.names = F,quote = F,sep = "\t")

head(df_str_length)
head(ref)
ref %>% select(STR,GC)
df_str_length %>% left_join(ref %>% select(STR,GC)) %>%
  ggplot(aes(x=TRGT,y=EH)) + 
  geom_point()

df_str_length %>% left_join(ref %>% select(STR,GC)) %>%
  ggscatter(.,x='TRGT',y="EH",color='GC',
          add = "reg.line",
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.y = 1000),
          xlab = "Long-read",
          ylab = "Short-read") + 
  geom_abline(slope = 1,linetype="dashed") + 
  theme_step1() + 
  theme(legend.position = 'right')
  

######
str_match_score_freqeuncy %>% filter(concordance == 0) %>% #filter(STR %in% dup_str$STR)
  count(MOTIFS) %>% #filter(n != 1) %>%
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(pct=n/sum(n)*100) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(MOTIFS = paste0(MOTIFS,'(',RU.length,')') ) %>%
  arrange(n, -RU.length) %>% #filter(n==2)
  mutate(MOTIFS = factor(MOTIFS, levels = unique(MOTIFS))) -> missmatch_count

#missmatch_count %>% write.table("figure2.allmissmatch_count.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  
head(missmatch_count)
missmatch_count %>% filter(n != 1) %>% select(MOTIFS)
missmatch_count[missmatch_count$n != 1,]$MOTIFS -> g1
missmatch_count[missmatch_count$n == 1,]$MOTIFS -> g2

g1
missmatch_count %>%
  ggplot(aes(x=MOTIFS,y=n,fill=factor(RU.length))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(n)),
            position=position_stack(vjust=0.5)) + 
  theme_step1() + 
  scale_x_discrete(limits = g1) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "None",
        plot.background = element_blank(),
        axis.text = element_text(size = 10)) + 
  guides(fill = guide_legend(reverse=F)) + 
  coord_flip(clip = 'off') #-> p1

p1  

missmatch_count %>%
  ggplot(aes(x=MOTIFS,y=n,fill=factor(RU.length))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(n)),
            position=position_stack(vjust=0.5)) + 
  theme_step1() + 
  scale_x_discrete(limits = g2) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "None",
        plot.background = element_blank(),
        axis.text = element_text(size = 10)) + 
  guides(fill = guide_legend(reverse=F)) + 
  scale_y_continuous(breaks = c(0, 1)) + 
  coord_flip(clip = 'off') #-> p1

missmatch_count %>%
  ggplot(aes(x=MOTIFS,y=n,fill=factor(RU.length))) + 
  geom_bar(stat='identity') + 
  theme_step1() + 
  guides(fill = guide_legend(reverse=F,title="RU length")) 

p3 <- get_legend(p3)


cowplot::plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(5,3,1)) -> p11
ggarrange(p2.1,p2.2,p3,nrow = 1, labels = c('b','c','d'), font.label = list(size = 28), label.y = 1.01,
          widths = c(1,1,2)) -> f1_up
gridExtra::grid.arrange(p11, bottom = "# of STRs",left = "Repeat Unit (RU length)")


### heatmap

library(ComplexHeatmap)
df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReaptnumber.with_dfiff.txt")

common_str_count <- read.table("figure1.frequency_common_motif.txt",header = T)
dup_str <- read.table("~/Desktop/KU/@research/STR/figure/rm_dup_str_list.txt",header = T)
head(common_str_count)
common_str_count %>% filter(n ==1) %>% dim()

head(df)
df %>% filter(!(STR %in% dup_str$STR)) %>% filter(STR %in% missmatch_str$STR) -> df_notmatch

df_notmatch %>% mutate(score = ifelse((STR1 == 0 & STR2 == 0),1,ifelse((STR1 != 0 & STR2 != 0),0,0.5))) %>% #head()
  select(ID,STR,MOTIFS,score) %>% group_by(ID,MOTIFS) %>%
  summarise(concordance = mean(score)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
  mutate(RU.length = str_length(MOTIFS)) -> forHeat

head(forHeat)
forHeat %>% filter(concordance == 0)

head(forHeat_toy)
forHeat %>% filter(MOTIFS %in% (common_str_count[common_str_count$n !=1,]$MOTIFS)) -> forHeat_toy
forHeat_toy %>% select(ID,MOTIFS,concordance) %>% pivot_wider(names_from = MOTIFS,values_from = concordance) -> forHeat_toy
forHeat_toy_m <- as.matrix(forHeat_toy)
head(forHeat_toy_m)
forHeat_toy_m[1:5,1:5]
forHeat_toy_m <-apply(forHeat_toy_m, 2, as.numeric)
rownames(forHeat_toy_m) <- forHeat_toy$ID
forHeat_toy_m[,2:ncol(forHeat_toy_m)] -> forHeat_toy_m
str_length(colnames(forHeat_toy_m))
str_length(gsub("[^CG]", "",colnames(forHeat_toy_m)))/str_length(colnames(forHeat_toy_m))
# Create annotations for RU.length and GC content
col_anno <- columnAnnotation(RU.length = str_length(colnames(forHeat_toy_m)), GC = str_length(gsub("[^CG]", "",colnames(forHeat_toy_m)))/str_length(colnames(forHeat_toy_m)))
ncol(forHeat_toy)
# Optionally, you can customize the annotation colors if needed
head(forHeat_toy_m)

# Generate the heatmap with clustering
Heatmap(forHeat_toy_m, 
        name = "Concordance Rate", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_dend = TRUE, 
        show_column_dend = TRUE, 
        row_dend_reorder = TRUE,
        show_column_names = FALSE,  # Hide x-axis text (sample names)
        show_row_names = FALSE,  
        column_dend_reorder = TRUE,
        top_annotation = col_anno)

