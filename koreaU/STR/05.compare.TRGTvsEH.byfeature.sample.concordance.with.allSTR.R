### compore length EH vs TRGT
library(tidyverse)
library(data.table)
library(ggbreak)
library(ggpubr)
library(ggtext)


df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReaptnumber.with_dfiff.txt")
str_match_score_freqeuncy <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt")
str_RU_count<- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.RUcount.txt",header =T)
str_qc_pass_common <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.txt",header = T)
ref <- read.table("/Users/ksmpooh/Desktop/KU/@research/STR/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed")
head(ref)
ref %>% mutate(ID=str_replace(str_split_fixed(V4,";",2)[,1],"ID=","")) %>% #head()
  select(-V4) ->ref
colnames(ref) <-c("chrom","start","end","STR")

head(str_RU_count)
head(str_qc_pass_common)
head(str_match_score_freqeuncy)
head(df)
df %>% count(ID) %>% count(n)
summary(df$nRU_diff_mean)


head(df)
head(df)
df %>% mutate(match = ifelse((STR1 == 0 & STR2 == 0),"ALL match",ifelse((STR1 != 0 & STR2 != 0),"ALL missmatch","One and One"))) %>% #-> df_match
  count(ID,match)-> df_match

df %>% mutate(match_bystr = ifelse((STR1 == 0 & STR2 == 0),"ALL match",ifelse((STR1 != 0 & STR2 != 0),"ALL missmatch","One and One"))) %>% #head()#-> df_match
  count(STR,MOTIFS,match_bystr) -> df_match_bystr

df %>% left_join(str_match_score_freqeuncy) %>% 
  mutate(EH = (EH_STR1+EH_STR2)/2) %>%
  mutate(TRGT = (TRGT_STR1+TRGT_STR2)/2) %>% select(-EH_STR1,-EH_STR2,-TRGT_STR1,-TRGT_STR2,-score) %>% #head()
  group_by(STR,MOTIFS) %>%
  summarise(EH = mean(EH),TRGT=mean(TRGT),concordance = mean(concordance),nRU_diff_mean= mean(nRU_diff_mean)) -> df_diff_meanMC_concordance



head(df_match)
head(df_match_bystr)

df_match %>% group_by(ID) %>% mutate(pct=n/sum(n)*100) 



df_match %>% group_by(ID) %>% mutate(pct=n/sum(n)*100) %>% #head()
  filter(match == 'ALL match') %>% #summary()
  filter(pct < 90) #%>% #summary()
  

df_match %>% group_by(ID) %>% mutate(pct=n/sum(n)*100)  %>% #head()
  #  mutate(n=factor(n)) %>%
  ggplot(aes(x=fct_reorder(ID,n),y=n,fill=fct_reorder(match,n))) +
  geom_bar(stat = 'identity') + 
  xlab("Sample") +
  ylab("# of STRs") + 
  scale_y_continuous(
    name = "# of STRs",
    sec.axis = sec_axis(~ . / max(302660) * 100, name = "Percentage of STR")
  ) +
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

head(df_diff_meanMC_concordance)
head(ref)

####  concordacne by chromosome
library(tidyverse)
library(ggtext)

head(ref)
head(str_match_score_freqeuncy)
head(ref)

chrOrder<-c(1:22,"X")
head(chrOrder)

head(df_diff_meanMC_concordance)

df_diff_meanMC_concordance %>% left_join(ref) %>% #head()
  mutate(chrom = str_replace(chrom,"chr","")) %>% #head()
  select(STR,MOTIFS,concordance,chrom,start) %>%
  mutate(chrom = factor(chrom, levels = chrOrder)) %>% 
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(pos = rank(start)) %>%
  ungroup() -> df_forMan

head(chrOrder)
head(df_forMan)

# 각 chromosome의 중앙 위치 계산 (axis_set)
axis_set <- df_forMan %>%
  group_by(chrom) %>%
  summarize(center = mean(pos))

head(df_forMan)

df_forMan %>% mutate(new_rankPOS = seq(1,nrow(df_forMan))) %>% 
  group_by(chrom) %>%
  summarize(center = mean(new_rankPOS)) -> axis_set

df_forMan %>% mutate(new_rankPOS = seq(1,nrow(df_forMan))) %>% #head()
  ggplot(aes(x = new_rankPOS, y = concordance, color = factor(chrom, levels = chrOrder))) +
  geom_point(alpha = 0.5,size=1) +
  scale_x_continuous(labels = axis_set$chrom,breaks = axis_set$center) +
  scale_color_manual(values = rep(
    c("#276FBF", "#183059"),
    length(unique(df_forMan$chrom))
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) + 
  labs(title = "STR Concordance by chromosome")  -> p1

p1
#  labs(x = "Chromosome", y = "Concordance",
#       title = "Manhattan Plot of STR Concordance (TRGT vs EH)") 


df_forMan %>% mutate(new_rankPOS = seq(1,nrow(df_forMan))) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
  ggplot(aes(x = new_rankPOS, y = concordance, color = GC)) +
  geom_point(alpha = 0.5,size=1) +
  scale_x_continuous(labels = axis_set$chrom,breaks = axis_set$center) +
  gradient_color(c("white","red")) + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) + 
  labs(title = "STR Concordance by GC (%)")  -> p2


df_forMan %>% mutate(new_rankPOS = seq(1,nrow(df_forMan))) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
  ggplot(aes(x = new_rankPOS, y = concordance, color = GC)) +
  geom_point(alpha = 0.5,size=1) +
  scale_x_continuous(labels = axis_set$chrom,breaks = axis_set$center) +
  gradient_color(c("white","red")) + 
  theme_minimal() +
  theme(
    legend.position = "bottom") -> p2_legend

p2_legend <- get_legend(p2_legend)


df_forMan %>% mutate(new_rankPOS = seq(1,nrow(df_forMan))) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>%
  ggplot(aes(x = new_rankPOS, y = concordance, color = RU.length)) +
  geom_point(alpha = 0.5,size=1) +
  scale_x_continuous(labels = axis_set$chrom,breaks = axis_set$center) +
  gradient_color(c("white","blue")) + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, vjust = 0.5)
  ) + 
  labs(title = "STR Concordance by RU length")  -> p3


df_forMan %>% mutate(new_rankPOS = seq(1,nrow(df_forMan))) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>%
  ggplot(aes(x = new_rankPOS, y = concordance, color = RU.length)) +
  geom_point(alpha = 0.5,size=1) +
  gradient_color(c("white","blue")) + 
  theme_minimal() +
  theme(legend.position = "bottom") + 
  labs(title = "STR Concordance by RU length")  -> p3_legend

p3_legend <- get_legend(p3_legend)


p4<-cowplot::plot_grid(p2_legend,p3_legend,nrow = 1)

p4

cowplot::plot_grid(p1,p2,p3,p4,nrow = 4,rel_heights = c(2,2,2,0.5)) + 
  labs(x = "Chromosome", y = "Concordance")


head(df_forMan)
chrOrder<-c(1:22,"X")
rev(chrOrder)

df_forMan %>% mutate(match = ifelse(concordance == 1,"ALL match",
                                    ifelse(concordance > 0.5, "0.5<x<1",
                                           ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  group_by(chrom,match) %>%
  count(chrom) %>% ungroup() %>% summarise(max = max(n)) -> max_n

df_forMan %>% mutate(match = ifelse(concordance == 1,"ALL match",
                                    ifelse(concordance > 0.5, "0.5<x<1",
                                           ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  group_by(chrom) %>%
  count(match) -> for_line

head(for_line)

df_forMan %>% mutate(match = ifelse(concordance == 1,"ALL match",
                                       ifelse(concordance > 0.5, "0.5<x<1",
                                              ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  #mutate(match = factor(match,levels=c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch"))) %>%
  mutate(match = factor(match, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  group_by(chrom) %>%
  count(match) %>% 
  ggplot(aes(x=factor(chrom,levels=rev(chrOrder)),y=n,fill=match)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = ifelse(match %in% c("ALL match"),n,""), y = 0),
            hjust = 0) +
  geom_text(aes(label = ifelse(match %in% c("0.5<x<1"),n,""),y = 0.75),
            hjust = 0) + 
  labs(x="Chromosome",y="Percentage") + 
  theme(legend.position= "none",
        legend.title = element_blank(),
   #     axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
  #      axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p1
  
df_forMan %>% mutate(match = ifelse(concordance == 1,"ALL match",
                                    ifelse(concordance > 0.5, "0.5<x<1",
                                           ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  #mutate(match = factor(match,levels=c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch"))) %>%
  mutate(match = factor(match, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  group_by(chrom,match) %>%
  count(chrom) %>%
  ggplot(aes(x=factor(chrom,levels=rev(chrOrder)),y=n,fill=match)) + 
  geom_bar(stat='identity') + 
  #geom_bar(stat='identity',position = 'fill') + 
  #geom_text(aes(label = ifelse(name %in% c("0.5<x<1","ALL match"),value,"")),
  geom_text(aes(label = ifelse(match %in% c("ALL match"),n,""), y = 0),
            hjust = 0) +
  geom_text(aes(label = ifelse(match %in% c("0.5<x<1"),n,"")),
            position = position_stack(vjust = 0.8)) + 
  labs(x="Chromosome",y="# of STRs") + 
  theme(legend.position = "none",
        legend.title = element_blank(),
             axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
              #axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p2

p2

df_forMan %>% mutate(match = ifelse(concordance == 1,"ALL match",
                                    ifelse(concordance > 0.5, "0.5<x<1",
                                           ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  mutate(match = factor(match, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  group_by(chrom,match) %>%
  count(chrom) %>%
  ggplot(aes(x=factor(chrom,levels=rev(chrOrder)),y=n,fill=match)) + 
  geom_bar(stat='identity') + 
  theme(legend.position = "bottom",
        legend.title = element_blank()) + 
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p3

p3 <- get_legend(p3)

p1
p2

cowplot::plot_grid(p1,p2,nrow = 1) -> p11
p11
#gridExtra::grid.arrange(p11,bottom = "Percentage",top="# of STRs",left = "Repeat Unit") -> p11

cowplot::plot_grid(p11,p3,nrow = 2,rel_heights = c(1,0.1))





df_forMan %>% mutate(match = ifelse(concordance == 1,"ALL match",
                                    ifelse(concordance > 0.5, "0.5<x<1",
                                           ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  #mutate(match = factor(match,levels=c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch"))) %>%
  mutate(match = factor(match, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  group_by(chrom) %>%
  count(match) %>% 
  mutate(cumulative_n = cumsum(n)) %>% #head()
  ggplot(aes(x=factor(chrom,levels=rev(chrOrder)),y=n,fill=match)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = ifelse(match %in% c("ALL match"),n,""), y = 0),
            hjust = 0) +
  geom_text(aes(label = ifelse(match %in% c("0.5<x<1"),n,""),y = 0.75),
            hjust = 0) + 
  labs(x="Chromosome",y="Percentage") + 
  theme(legend.position= "none",
        legend.title = element_blank(),
        #     axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #      axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() +
  geom_line(aes(y = n/max_n$max, group = match), color = match, size = 1) #+
  scale_y_continuous(sec.axis = sec_axis(~ .*n/max_n$max, name = "Cumulative Count")) 


#######
ggplot


  geom_text(aes(label = ifelse(name %in% c("0.5<x<1"),value,""), y = 0.5),
            hjust = 0) + 
  geom_text(aes(label = ifelse(name %in% c("ALL missmatch"),value,""),fontface = "bold",y = 1), 
            hjust = 0,color = "darkred",) + 
  geom_text(aes(label = ifelse(name %in% c("0<x=<0.5"),add_text,"")), 
            position = position_fill(vjust = 0),color='darkgreen') + 
  scale_y_continuous(labels = percent) + 
  #labs(title = "Other, Rank 1~19") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip()






#########


head(df)
head(str_match_score_freqeuncy)

head(df_diff_meanMC_concordance)
head(str_RU_count)
colnames(str_RU_count)[1] <- "MOTIFS"

df_diff_meanMC_concordance %>% #group_by(MOTIFS) %>% #summarise(EH = mean(EH),TRGT = mean(TRGT)) %>%  
  left_join(str_RU_count) %>% #head()
  mutate(count=ifelse(n==1,"n=1",
                      ifelse(n==2,"n=2",
                             ifelse(Rank %in% c(1:10),"Top10","other")))) %>% #ungroup() %>% count(count)
  #filter(count == "Top10")
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
    mutate(RU.length = str_length(MOTIFS)) %>% #head()
  ggscatter(.,x='TRGT',y="EH",color = "GC",size='RU.length',
            #  ggscatter(.,x='mean_MC_TRGT',y="mean_MC_EH",size = "GC",color='RU.length',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            #xlim = c(0,0.5),
            #ylim = c(0,0.5),
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n",label.y.npc = "top"),
                        #cor.coeff.args = list(method = "pearson", label.sep = "\n"),
            xlab = "Long-read (TRGT)",
            ylab = "Short-read (EH)") + 
  geom_abline(slope = 1,linetype="dashed") + 
  gradient_color(c("blue","red")) + 
  #  guides(color=guide_legend(title="GC contents of RU")) + 
  theme(legend.position = "right",strip.text = element_text(size = 20)) + 
  facet_wrap(~count,nrow = 2,scales = 'free') -> p1

p1

df_diff_meanMC_concordance %>% #head()
  mutate(concordance_group = ifelse(concordance == 1,"ALL match",
                                    ifelse(concordance > 0.5, "0.5<x<1",
                                           ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>% filter(concordance == 1)

df_diff_meanMC_concordance %>% #head()
  mutate(concordance_group = ifelse(concordance == 1,"ALL match",
                              ifelse(concordance > 0.5, "0.5<x<1",
                                     ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  mutate(concordance_group = factor(concordance_group,levels=c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch"))) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
  mutate(RU.length = str_length(MOTIFS)) %>% #head()
  #mutate(concordance = as.factor(concordance)) %>%
  ggscatter(.,x='TRGT',y="EH",color = "GC",size='RU.length',
            #  ggscatter(.,x='mean_MC_TRGT',y="mean_MC_EH",size = "GC",color='RU.length',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            #xlim = c(0,0.5),
            #ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n"),
            xlab = "Long-read (TRGT)",
            ylab = "Short-read (EH)") + 
  geom_abline(slope = 1,linetype="dashed") + 
  #gradient_color(c("white","red")) + 
  gradient_color(c("blue","red")) + 
  #  guides(color=guide_legend(title="GC contents of RU")) + 
  theme(legend.position = "right",strip.text = element_text(size = 20)) + 
  facet_wrap(~concordance_group,scales = 'free') -> p2

  
cowplot::plot_grid(p1,p2)



## by length
ru.scale = c(seq(1,14),"15+")
ru.scale
ru %>% select(MOTIFS) %>% mutate(RU.length = str_length(MOTIFS)) %>%
  count(RU.length)



df_diff_meanMC_concordance %>% #group_by(MOTIFS) %>% #summarise(EH = mean(EH),TRGT = mean(TRGT)) %>%  
  left_join(str_RU_count) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
  mutate(RU.length = str_length(MOTIFS)) %>% #count(RU.length)
  ungroup() %>%
  mutate(RU.length = ifelse(RU.length >= 15,"15+",RU.length)) -> plot.data

cor_with_p <- function(x, y) {
  test <- cor.test(x, y, method = "pearson")
  tibble(correlation = test$estimate, p_value = test$p.value)
}
plot.data %>% filter(RU.length=="2") -> a
a
cor.test(a$EH,a$TRGT,method = "pearson")
cor(a$EH,a$TRGT,method = "pearson")

#  head()
plot.data %>% group_by(RU.length) %>% 
  summarise(R = cor(TRGT, EH, method = "pearson"),P=)

head(plot.data)
plot.data %>% group_by(RU.length) %>%
  summarize(result = list(cor_with_p(TRGT, EH))) %>%
  unnest(result) %>% mutate(p_value = ifelse(p_value == 0,2.2e-16,p_value)) -> plot.data.cor

head
plot.data %>% 
  ggscatter(.,x='TRGT',y="EH",color = "GC",#size='RU.length',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.sep = "\n",label.y=400),
            xlab = "Long-read (TRGT)",
            ylab = "Short-read (EH)") + 
  geom_abline(slope = 1,linetype="dashed") + 
  #stat_cor(method = "pearson", label.x.npc = "center", label.y.npc = "top") + 
  gradient_color(c("blue","red")) + 
  #  guides(color=guide_legend(title="GC contents of RU")) + 
  theme(legend.position = "right",strip.text = element_text(size = 10)) + 
  facet_wrap(~factor(RU.length,levels = ru.scale),nrow = 2) -> p1

p1

plot.data %>% #left_join(plot.data.cor)
  ggscatter(.,x='TRGT',y="EH",color = "GC",#size='RU.length',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.sep = "\n",
                                  label.x.npc = "center",label.y.npc = "top"),
            xlab = "Long-read (TRGT)",
            ylab = "Short-read (EH)") + 
  geom_abline(slope = 1,linetype="dashed") + 
  #stat_cor(method = "pearson", label.x.npc = "center", label.y.npc = "top") + 
  gradient_color(c("blue","red")) + 
  #  guides(color=guide_legend(title="GC contents of RU")) + 
  theme(legend.position = "right",strip.text = element_text(size = 10)) + 
  facet_wrap(~factor(RU.length,levels = ru.scale),nrow = 2,scale='free') -> p2


p1
p2


###################

head(df_diff_meanMC_concordance)



df %>% ggplot(aes(x=ID,y=nRU_diff_mean,fill=nRU_diff_mean)) + 
  geom_violin()+
  scale_y_break(c(50,300))+
  #scale_y_break(c(700, 4900),ticklabels = c(1,2,3,4950))
  theme(axis.text.x = element_blank())


head(max(df_match$n))
head(df_match)


df %>% mutate(match_count = ifelse((STR1 == 0 & STR2 == 0),2,ifelse((STR1 != 0 & STR2 != 0),0,1))) %>% #head()
  group_by(STR,MOTIFS) %>% summarise(match_count = sum(match_count)) -> df_match_count

head(df)
head(df_match_count)

df %>% filter(ID != "NIH23J3904558") %>%
  filter(STR1 == 0,STR2 == 0)  %>% select(ID,STR,MOTIFS) %>% count(STR,MOTIFS) -> a_65
df %>% filter(STR1 == 0,STR2 == 0)  %>% select(ID,STR,MOTIFS) %>% count(STR,MOTIFS) -> a

df %>% filter(ID != "NIH23J3904558") %>% filter(STR1 != 0 | STR2 != 0) %>% 
  count(STR,MOTIFS)  -> b_65
  
df %>% filter(STR1 != 0 | STR2 != 0) %>% 
  count(STR,MOTIFS)  -> b

df %>% filter(ID != "NIH23J3904558") %>% filter(STR1 != 0,STR2 != 0) %>% 
  count(STR,MOTIFS)  -> c_65

df %>% filter(STR1 != 0,STR2 != 0) %>% 
  count(STR,MOTIFS)  -> c


head(c)
head(a)
head(b)
a %>% filter(n==66) %>% dim
a_65 %>% filter(n==65) %>% dim

b %>% dim
b_65 %>% dim

c %>% filter(n==66) %>% dim
c_65 %>% filter(n==65) %>% dim


a %>% rbind(b) %>% rbind(c) %>% filter(n ==0)


#write.table(a,"~/Desktop/KU/@research/STR/02.compare/STR.longvsshort.match.count.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#a <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.longvsshort.match.count.txt")

head(df_match_count)
dim(df_match_count)
head(df)
#ifelse((STR1 == 0 & STR2 == 0),"ALL match",ifelse((STR1 != 0 & STR2 != 0),"ALL missmatch","One and One"))) %>% #-> df_match
#"ALL match","ALL missmatch","One and One"
df %>% mutate(match_bystr = ifelse((STR1 == 0 & STR2 == 0),"match",'missmatch')) %>% #head()
  group_by(STR,MOTIFS) %>%
  #summarise(sum = count(match_bystr))-> test
  count(STR,MOTIFS,match_bystr) #-> df_match_bystr



df_match_bystr %>% group_by(match_bystr) %>%
  count(n) %>%
  ggplot(aes(x=n,y=nn,fill=factor(match_bystr,levels = c("ALL missmatch","One and One","ALL match")))) +
  geom_bar(stat='identity') + 
  xlab("# of Sample") +
  ylab("# of STR") +
  scale_y_break(c(51000,140000),ticklabels = c(0,10000,20000,30000,40000,50000,145000))+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 10)) + 
  facet_wrap(~match_bystr,ncol = 3)
  #facet_wrap(~match_bystr,ncol = 3,scales = 'free_y')


df_match_bystr %>% group_by(match_bystr) %>% filter(n==66) %>%
  count(match_bystr)
  


df_match_bystr %>% mutate(score = ifelse(match_bystr == "ALL match",n*2,if_else(match_bystr=="ALL missmatch",n*0,n*1))) %>%
  group_by(STR,MOTIFS) %>% #head()
  summarise(score=sum(score)) %>%
  mutate(concordance = score/132) %>%
  ungroup() -> str_match_score_freqeuncy

head(str_match_score_freqeuncy)

#write.table(str_match_score_freqeuncy,"~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df_match_bystr %>% mutate(score = ifelse(match_bystr == "ALL match",n*2,if_else(match_bystr=="ALL missmatch",n*0,n*1))) %>%
  group_by(STR,MOTIFS) %>% #head()
  summarise(score=sum(score)) %>% #head()
  ungroup() %>%
  group_by(MOTIFS) %>%
  mutate(allelenumber=132) %>%
  summarise(score=sum(score),total_allele=sum(allelenumber)) %>%
  mutate(concordance = score/total_allele) %>% 
  ungroup() -> str_match_score_freqeuncy_byRU

head(str_match_score_freqeuncy_byRU)
str_match_score_freqeuncy_byRU %>% 
  


head(df)
df %>% select(ID,STR,EH_STR1,EH_STR2,TRGT_STR1,TRGT_STR2) %>% 
  mutate(EH_mean = (EH_STR1+EH_STR2)/2) %>%
  mutate(TRGT_mean = (TRGT_STR1+TRGT_STR2)/2) %>% select(EH_STR1,STR2)



  #count()
head(test)
df_match_bystr %>% group_by(match_bystr)%>%summarise(sum = sum(n))


df_match_bystr %>% filter(n == 66) %>% count(match_bystr)

df_match_bystr %>% count(match_bystr,n) %>% #head()
  ggplot(aes(x=n,y=nn,fill=factor(match_bystr,levels = c("ALL missmatch","One and One","ALL match")))) +
  geom_bar(position = 'stack',stat='identity') + 
  xlab("# of Sample") +
  ylab("# of STR") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 10))


df_match_bystr %>% count(match_bystr,n) %>% filter(nn == 0)
  
## str1 | str2 둘중 하나라도 0인거
b %>% count(n) %>% #head()
  ggplot(aes(x=n,y=nn)) + 
  geom_bar(stat='identity') + 
  scale_y_break(c(20000,45000),ticklabels = c(0,5000,10000,15000,20000,45000,50000))+
  xlab("# of Sample") +
  ylab("# of STR")

## str1, str2 둘다 0인거
c %>% count(n) %>% #head()
  ggplot(aes(x=n,y=nn)) + 
  geom_bar(stat='identity') + 
  scale_y_break(c(10000,30000),ticklabels = c(0,10000,30000,31000,32000,33000,34000))+
  xlab("# of Sample") +
  ylab("# of STR")

## b+c
b %>% count(n) %>% left_join(c %>% count(n) %>% rename(`All allele not match` = nn)) %>%
  mutate(`One of alleles not match` = nn-`All allele not match`) %>% #head()
  select(-nn) %>% pivot_longer(cols = 2:3) %>%
  ggplot(aes(x=n,y=value,fill=name)) + 
  geom_bar(position = 'stack',stat='identity') + 
  scale_y_break(c(20000,45000),ticklabels = c(0,5000,10000,15000,20000,45000,50000)) +
  xlab("# of Sample") +
  ylab("# of STR") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 10))

#ggplot(aes())

#### all match RU type count
c %>% filter(n==66) %>% dim()

head(a)

head(b)
head(a)
head(c)

a %>% filter(n == 66) %>% 
  count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) -> str_RU_count

a %>% filter(n == 66) -> all_sample_same_RUtype
head(str_RU_count)

str_RU_count %>% #head()
  mutate(Rank = ifelse(Rank %in% c(1:10),Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Rank %in% c(1:10),MOTIFS,"Other")) %>% 
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(pct=n/sum(n)*100) %>% #head()
  ggplot(aes(x=reorder(MOTIFS,-n),y=n,fill=reorder(MOTIFS,-n))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(round(pct,1), '%')),
            position=position_stack(vjust=0.99)) + 
  #  scale_y_cut(breaks=c(50000, 160000), scales=c(0)) +
  scale_y_break(c(15000, 94000),ticklabels = c(0,5000,10000,15000,94000,96000)) +
  ylab("# of STR") + 
  theme(axis.title.x = element_blank(),
        legend.position = "None",
        axis.text = element_text(size = 10)) -> p1
p1
head(c)

c %>% filter(n == 66) %>% #dim()
  count(MOTIFS) %>% arrange(-n) %>% #ggplot(aes(x=MOTIFS,y=n)) + geom_bar(stat = 'identity')
#  mutate(Rank = ifelse(n == 1,'Other',n)) %>% 
  mutate(MOTIFS = ifelse(n == 1,'n=1',MOTIFS)) %>%
  mutate(MOTIFS = ifelse(n == 2,'n=2',MOTIFS)) %>% 
  mutate(MOTIFS = ifelse(n == 3,'n=3',MOTIFS)) %>%
  mutate(MOTIFS = ifelse(n == 4,'n=4',MOTIFS)) %>% 
  mutate(MOTIFS = ifelse(n == 5,'n=5',MOTIFS)) %>% 
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(pct=n/sum(n)*100) %>% #head()
  ggplot(aes(x=reorder(MOTIFS,-n),y=n,fill=reorder(MOTIFS,-n))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(round(pct,1), '%')),
            position=position_stack(vjust=0.99)) + 
  ylab("# of STR") + 
  theme(axis.title.x = element_blank(),
        legend.position = "None",
        axis.text = element_text(size = 10)) -> p2

p2

b %>% filter(n == 66) %>% #dim()
  count(MOTIFS) %>% arrange(-n) %>% #ggplot(aes(x=MOTIFS,y=n)) + geom_bar(stat = 'identity')
  mutate(Rank = rank(-n)) %>% 
  mutate(Rank = ifelse(Rank %in% c(1:10),Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Rank %in% c(1:10),MOTIFS,"Other")) %>% 
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(pct=n/sum(n)*100) %>% #head()
  ggplot(aes(x=reorder(MOTIFS,-n),y=n,fill=reorder(MOTIFS,-n))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(round(pct,1), '%')),
            position=position_stack(vjust=0.99)) + 
  ylab("# of STR") + 
  theme(axis.title.x = element_blank(),
        legend.position = "None",
        axis.text = element_text(size = 10))

head(df)
head(df_match_bystr)
head(str_RU_count)


df %>% filter(ID == "NIH23F1724125") %>% select(STR,MOTIFS) %>% unique() %>% #head()
  count(MOTIFS) %>% arrange(-n) %>% mutate(Ori_Rank = rank(-n)) -> str_all_count

head(str_all_count)
head(str_RU_count)
sum(str_RU_count$n)
colnames(str_RU_count)[2] <- "match_n"
colnames(str_RU_count)[3] <- "match_Rank"

str_all_count %>% left_join(str_RU_count) %>% mutate(notmatch=n-match_n) %>% 
  select(MOTIFS,Ori_Rank,match_Rank,match_n,notmatch) %>% #head()
  #select(-match_Rank) %>%
  mutate(match_Rank = ifelse(Ori_Rank %in% c(1:10),match_Rank,"Other")) %>%
  mutate(Ori_Rank = ifelse(Ori_Rank %in% c(1:10),Ori_Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Ori_Rank %in% c(1:10),MOTIFS,"Other")) %>% na.omit() %>% #head()
  group_by(MOTIFS,Ori_Rank,match_Rank) %>% summarise(match_n = sum(match_n),notmatch=sum(notmatch)) %>%
  pivot_longer(cols = c(match_n,notmatch)) %>% #count(MOTIFS)
  mutate(name = ifelse(name == "match_n","ALL match","missmatch")) %>% #head()
  mutate(name = factor(name,levels = c("ALL match","missmatch"))) %>%  #head
  group_by(MOTIFS) %>% mutate(pct=value/sum(value)*100) -> plot_data

plot_data %>%
  ggplot(aes(x=reorder(MOTIFS,-value),y=value,fill=name)) + 
  geom_bar(stat = 'identity',position = 'stack') + 
  ylab("# of STR") + 
  labs(title = "Percentage of Motif Counts for ALL Matches",subtitle = "(RU count) Original Rank -> ALL match Rank") + 
  geom_text(aes(label = ifelse(name == "ALL match", paste0(round(pct, 1), '%'), '')),
            position = position_stack(vjust = 0.99)) +
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,size=20),
        plot.subtitle = element_text(hjust = 0.5,size=10),
        legend.title = element_blank(),
        axis.text = element_text(size = 10)) -> p1

p1 + annotate("text", x = plot_data$MOTIFS, y = -10, # y 값을 적절히 조정
             label = ifelse(plot_data$name == "missmatch" & plot_data$Ori_Rank != "Other", paste0(plot_data$Ori_Rank, '->', plot_data$match_Rank), '')) -> p1
   
p1
head(plot_data)


geom_bar(position="fill", stat="identity")


plot_data %>% filter(match_Rank != "Ohter") %>% select(MOTIFS,Ori_Rank,match_Rank)

a %>% filter(n == 66) %>% 
  count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) -> str_RU_count

a %>% filter(n == 66) -> all_sample_same_RUtype
head(str_RU_count)

str_RU_count %>% #head()
  mutate(Rank = ifelse(Rank %in% c(1:10),Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Rank %in% c(1:10),MOTIFS,"Other")) %>% 
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(pct=n/sum(n)*100) %>% #head()
  ggplot(aes(x=reorder(MOTIFS,-n),y=n,fill=reorder(MOTIFS,-n))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(round(pct,1), '%')),
            position=position_stack(vjust=0)) + 
  #  scale_y_cut(breaks=c(50000, 160000), scales=c(0)) +
  scale_y_break(c(15000, 94000),ticklabels = c(0,5000,10000,15000,94000,96000)) +
  ylab("# of STR") + 
  theme(axis.title.x = element_blank(),
        legend.position = "None",
        axis.text = element_text(size = 10))

head(c)





p1
p2
#cowplot::plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(1,1,0.2))

head(a)
head(b)
head(str_RU_count)
head(all_sample_same_RUtype)

##########
head(df)
df %>% mutate(EH = (EH_STR1 + EH_STR2)/2,TRGT = (TRGT_STR1 + TRGT_STR2)/2) %>%
  group_by(STR,MOTIFS) %>% summarise(mean_EH=mean(EH),mean_TRGT=mean(TRGT)) -> df_str_mean_repeat_count

head(str_match_score_freqeuncy)
str_match_score_freqeuncy %>% 
  mutate(RU.length = str_length(MOTIFS)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>% 
  ggplot(aes(x=as.factor(RU.length),y=concordance)) + 
  xlab("RU.length") + 
  ylab("concordance by STR") + 
  geom_violin() -> p1
p1

str_match_score_freqeuncy_byRU %>% #head()
  #filter(concordance != 1) %>%
  mutate(RU.length = str_length(MOTIFS)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>% 
  ggplot(aes(x=as.factor(RU.length),y=concordance)) + 
  geom_violin() +
  xlab("RU.length") + 
  ylab("concordance by RU")-> p2

cowplot::plot_grid(p1,p2,nrow = 2)


str_match_score_freqeuncy %>% mutate(g = ifelse(concordance==1,"ALL match","missmatch")) %>%
  mutate(g = ifelse(concordance==0,"ALL missmatch",g)) %>% #head()
  count(g)

str_match_score_freqeuncy %>% mutate(RU.length = str_length(MOTIFS)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>% #head()
  mutate(g = ifelse(concordance==1,"ALL match (148,514)","missmatch (153,963)")) %>%
  mutate(g = ifelse(concordance==0,"ALL missmatch (189)",g)) %>% #head()
  select(STR,RU.length,g) %>% 
  ggplot(aes(x=RU.length,fill=g)) + 
  geom_density(alpha=0.7) +
  theme(legend.position = 'none')-> p1
  
str_match_score_freqeuncy %>% mutate(RU.length = str_length(MOTIFS)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>% #head()
  mutate(g = ifelse(concordance==1,"ALL match (148,514)","missmatch (153,963)")) %>%
  mutate(g = ifelse(concordance==0,"ALL missmatch (189)",g)) %>% #head()
  select(STR,GC,g) %>% 
  ggplot(aes(x=GC,fill=g)) + 
  geom_density(alpha=0.7) +
  labs(fill="STR match type(# of STR)") -> p2
  
cowplot::plot_grid(p1,p2)

####### total STR
################################################## check_~20240710
head(df)
df %>% mutate(match_type = ifelse(STR %in% b$STR,"missmatch","match")) %>% 
  mutate(match_type=ifelse(STR %in% c$STR,"All_missmatch",match_type)) -> df_type

df_type %>% mutate(STR1_length = TRGT_STR1*str_length(MOTIFS),STR2_length = TRGT_STR2*str_length(MOTIFS)) %>%
  #filter(STR %in% c1$STR) %>% 
  mutate(RU_length = str_length(MOTIFS)) %>%
  mutate(RU_length = ifelse(RU_length %in% c(1:6),RU_length,"7~"))-> df_type

## 이거 오래 걸림
head(df_type)
df_type %>% #filter(ID == "NIH23F1724125" )%>% 
  ggpubr::ggscatter(.,x='STR1_length',y="STR2_length",#shape = 21,
                    #fill = "match_type",
                    color = "match_type",#,color='Rank',
                    add = "reg.line",
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #xlim = c(0,0.5),
                    #ylim = c(0,0.5),
                    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                    #cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                    xlab = "STR1_length",
                    ylab = "STR2_length") + 
  geom_abline(slope = 1,linetype="dashed") + 
  theme(legend.position = "bottom") + 
  facet_wrap(~RU_length,nrow = 2)


#######
#head(new_df)
######### non-match


head(df)
df %>% filter(STR %in% all_sample_same_RUtype$STR) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(mean_repeat_count = (TRGT_STR1 + TRGT_STR2)/2) %>%
  mutate(Whole.length = RU.length*mean_repeat_count) %>% #head()
  group_by(STR,MOTIFS) %>%
  summarise(RU.length = mean(RU.length),mean_Whole.length = mean(Whole.length),mean_repeat_count= mean(mean_repeat_count)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> new_df

head(str_match_score_freqeuncy)
head(new_df)
new_df %>% left_join(str_match_score_freqeuncy) %>% #head()
  ggplot(aes(x=RU.length,y=mean_repeat_count,color=GC)) +
  geom_point() +
  #theme_light() + 
  scale_color_gradient(low = "white", high = "red")




head(a)
a %>% count(n)
#str_RU_count
head(all_sample_same_RUtype)
dim(all_sample_same_RUtype)



new_df %>% 
  #head()
  
  ggplot(aes(x=RU.length,y=mean_repeat_count,color=GC,fill=)) +
  geom_point() +
  #theme_light() + 
  scale_color_gradient(low = "white", high = "red")





#head(all_sample_same_RUtype)
df %>% filter(STR %in% all_sample_same_RUtype$STR) %>% select(ID,STR,MOTIFS,TRGT_STR1,TRGT_STR2) %>% 
  mutate(compare_repeat_count_A1vsA2 = ifelse(TRGT_STR1 == TRGT_STR2,"A1==A2","A1 != A2")) %>% #head
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(mean_repeat_count = (TRGT_STR1 + TRGT_STR2)/2) %>%
  group_by(STR,MOTIFS,compare_repeat_count_A1vsA2) %>%
  summarise(RU.length = mean(RU.length),mean_repeat_count= mean(mean_repeat_count)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>% #head()
  ungroup() %>%
  ggplot(aes(x=as.factor(RU.length),y=mean_repeat_count,fill=compare_repeat_count_A1vsA2)) +
  geom_boxplot() +
  xlab("RU.length") +
  ylab("Mean of repeat count") +
  guides(fill=guide_legend(title="Repeat count comparison")) + 
  #guide_legend(title="Repeat Count compare (A1 == A2)") + 
  theme(legend.position = "bottom")
#facet_grid(~compare_repeat_count_A1vsA2)



str_class_ref <- read_table("~/Desktop/KU/@research/STR/db/hg38_2020_rmsk.bed",col_names = F)
#ref <- read_table("~/Desktop/KU/@research/STR/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.rmdup.bed",col_names = F)
head(ref)
head(str_class_ref)

ref %>% select(X1,X2,X3,X6,X8) -> ref
str_class_ref %>% mutate(ID = paste0(X1,"_",X2,"_",X3)) -> str_class_ref

head(ref)
head(str_class_ref)

str_class_ref %>% filter(ID %in% ref$X8) %>% count(X4) -> x

head(df)
df %>% select(STR) %>% unique() -> x1
head(x)
head(x1)

table(x1$STR %in% str_class_ref$ID)
head(str_class_ref_common)

all_sample_same_RUtype

head(str_class_ref)
head(all_sample_same_RUtype)
head(df)
##
# eh
#load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.info.RData")
load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.PASS.intersect.RData")

#head(eh_pass_common)
head(eh)
eh %>% filter(FILTER == "PASS") %>% filter(STR_ID %in% df$STR) -> eh

head(b)
head(c)
head(eh)
eh %>% mutate(match_type = ifelse(STR_ID %in% b$STR,"missmatch","match")) %>% 
  mutate(match_type=ifelse(STR_ID %in% c$STR,"All_missmatch",match_type)) -> eh_pass_common_matchtype_QCinfo


head(eh_pass_common_matchtype_QCinfo)
eh_pass_common_matchtype_QCinfo %>% select(-FILTER,-type) %>% 
  mutate(ADSP1 = as.integer(str_split_fixed(ADSP,"/",2)[,1]),ADSP2 = as.integer(str_split_fixed(ADSP,"/",2)[,2])) -> eh_pass_common_matchtype_QCinfo_adsp

head(eh_pass_common_matchtype_QCinfo_adsp)

rm(eh)
rm(eh_pass_common_matchtype_QCinfo)
rm(eh_pass_common_matchtype_QCinfo_adsp)
#rm(str_)


head(df)
df %>% mutate(match_type = ifelse(STR %in% b$STR,"missmatch","match")) %>% 
  mutate(match_type=ifelse(STR %in% c$STR,"All_missmatch",match_type)) -> df_type

head(df_type)
df_type %>% filter(!match_type %in% c("match"))
head(c);dim(c)

c %>% filter(n == 66) -> c1
head(c)
head(c1)
head(df_type)
df_type %>% mutate(STR1_length = TRGT_STR1*str_length(MOTIFS),STR2_length = TRGT_STR2*str_length(MOTIFS)) %>%
  #filter(STR %in% c1$STR) %>% 
  mutate(RU_length = str_length(MOTIFS)) %>%
  mutate(RU_length = ifelse(RU_length %in% c(1:6),RU_length,"7~"))-> df_type


head(df_type)
df_type %>% #filter(ID == "NIH23F1724125" )%>% 
  ggpubr::ggscatter(.,x='STR1_length',y="STR2_length",#shape = 21,
                    #fill = "match_type",
                    color = "match_type",#,color='Rank',
                    add = "reg.line",
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #xlim = c(0,0.5),
                    #ylim = c(0,0.5),
                    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                    #cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                    xlab = "STR1_length",
                    ylab = "STR2_length") + 
  geom_abline(slope = 1,linetype="dashed") + 
  theme(legend.position = "bottom") + 
  facet_wrap(~RU_length,nrow = 2)
#facet(RU_length,ncol = 3)
#facet_grid(~RU_length,)



ggplot(aes(x=STR1_length,y=STR2_length,color=RU_length)) +
  geom_point() + 
  facet_grid(~match_type)

head(ref)

ref %>% select(6:17)
head(df)
df %>% select(STR, MOTIFS) %>% unique() -> str_type
head(str_type)
head(ref)

ref %>% select(6,7,8,9,10,11,12,13) -> ref
ref %>% mutate(STR = paste0(X6,"_",X7,"_",X8)) -> ref
head(ref)
table(str_type$STR %in% ref$X11)
str_type %>% left_join(ref) -> str_type
#common_str <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/TRID_common_EHpass_TRGTupper0.8.txt",header = T)

#str_type %>% mutate(chrom = str_split_fixed(STR,"_",3)[,1],start = str_split_fixed(STR,"_",3)[,2],end = str_split_fixed(STR,"_",3)[,3]) ->str_type
str_type %>% select(STR,MOTIFS,chrom,start,end) ->str_type
head(str_type)
#str_type %>% filter(start != "") -> str_type
#write.table(str_type,"/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.txt",col.names = T,row.names = F,quote = F,sep = "\t")

head(str_type)
str_type %>% filter(is.na(X10)) %>% count()
str_type %>% na.omit() %>% select(STR,MOTIFS,X12)->str_type
heads(str_type)
head(df)
head(df_type)

df_type %>% left_join(str_type) %>% na.omit() -> str_type
head(str_type$ID)
table(str_type$X12)
str_type %>% filter(ID == "NIH23F1724125") %>%
  select(STR,MOTIFS,match_type,RU_length,X12) %>% unique() -> str_type1

head(str_type1)

str_type1 %>% count(match_type,X12) %>%
  pivot_wider(names_from = match_type,values_from = n)


#count(match_type,X12,)
str_type %>% 
  ggplot(aes(x=))



#####


df_match_count %>% ungroup() %>% #group_by(STR,MOTIFS) %>%
  mutate(match_count_freq = match_count/max(match_count)) %>%  
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
  mutate(RU.length = str_length(MOTIFS)) %>% #head()
  mutate(RU.length = ifelse(RU.length %in% c(1:10),RU.length,"10~")) %>% #head()
  ggplot(aes(x=RU.length,y=match_count_freq,fill=GC)) +
  geom_violin()
###

#load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.upper0.8.intersect.RData")
#load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.info.RData")


head(df_match_count)
head(df)
head(tgrt)
tgrt %>% filter(TRID %in% df_match_count$STR) %>% 
  select(-STRUC,-MC,-AP) %>% pivot_wider(names_from = Allele,values_from = AM) %>% 
  filter(allele_1 != "." | allele_2 != ".") -> tgrt
colnames(tgrt) <- c()



trgt %>% group_by(TRID)
#####