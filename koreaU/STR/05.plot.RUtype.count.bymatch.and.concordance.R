###### compare long vs short STR type motif compare

library(tidyverse)
library(data.table)
library(ggbreak)
library(ggpubr)
library(ggtext)
library(scales)

str_match_score_freqeuncy <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt")
str_RU_count<- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.RUcount.txt",header =T)
str_qc_pass_common <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.txt",header = T)
head(str_RU_count)
colnames(str_RU_count)[1] <- "MOTIFS"

#head(str_match_score_freqeuncy)
str_match_score_freqeuncy %>% mutate(match = ifelse(concordance == 1,"Match","Missmatch")) %>% 
  group_by(MOTIFS,match) %>%
  count(MOTIFS) %>% #head()
  pivot_wider(names_from = match,values_from = n) -> str_RU_count_bymatch
  

head(str_RU_count)
head(str_RU_count_bymatch)

head(str_RU_count)

rank_min = 1
rank_max = 99

str_RU_count %>% left_join(str_RU_count_bymatch) %>% #head()
  select(-n) %>% #head()
  pivot_longer(cols = c("Match","Missmatch")) %>% #head()
  mutate(Rank = ifelse(Rank %in% c(rank_min:rank_max),Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Rank %in% c(rank_min:rank_max),MOTIFS,"Other")) %>% na.omit() %>% #head()
  group_by(MOTIFS,Rank,name) %>% 
  summarise(value = sum(value)) %>% 
  left_join(str_RU_count %>% select(-Rank)) %>% #head()
  filter(MOTIFS != "Other") %>% mutate(RU.length = str_length(MOTIFS))-> plot.data

str_RU_count %>% left_join(str_RU_count_bymatch) %>% #head()
    pivot_longer(cols = c("Match","Missmatch")) %>% #head()
  mutate(Rank = ifelse(Rank %in% c(rank_min:rank_max),Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Rank %in% c(rank_min:rank_max),MOTIFS,"Other")) %>% na.omit() %>% 
  filter(Rank == "Other") %>%
  group_by(MOTIFS,Rank,name) %>%
  summarise(value = sum(value)) %>% ungroup() %>%
  mutate(n = sum(value)) -> plot.data.other

head(plot.data)
head(plot.data.other)

head(plot.data %>% filter(Rank %in% c(1:19)))
plot.data %>% filter(Rank %in% c(1:19)) %>% rbind(plot.data.other)%>% #filter(MOTIFS == "Other")
  mutate(name = factor(name, levels = c("Missmatch", "Match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = value), 
            position = position_fill(vjust = 0.5)) + 
  labs(title = "Other, Rank 1~19") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        #axis.text.x = element_text(colour = RU.length),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p1

p1

head(rank_range)
rank_range

rank_range <- as.character(20:39)
rank_pattern <- paste(rank_range, collapse = "|")
rank_pattern


plot.data %>% 
  filter(str_detect(Rank, rank_pattern)) %>% #head
  mutate(name = factor(name, levels = c("Missmatch", "Match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = value), 
            position = position_fill(vjust = 0.5)) + 
   labs(title = "Rank 20~39") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p2

p2
p1

rank_range <- as.character(40:59)
rank_pattern <- paste(rank_range, collapse = "|")
rank_pattern

plot.data %>% 
  filter(str_detect(Rank, rank_pattern)) %>% #head
  mutate(name = factor(name, levels = c("Missmatch", "Match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = value), 
            position = position_fill(vjust = 0.5)) + 
  labs(title = "Rank 40~59") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) +
  coord_flip() -> p3

p2
p3

rank_range = c(60:79)
rank_pattern <- paste(rank_range, collapse = "|")
rank_pattern

plot.data %>% 
  filter(str_detect(Rank, rank_pattern)) %>% #head
  mutate(name = factor(name, levels = c("Missmatch", "Match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = value), 
            position = position_fill(vjust = 0.5)) + 
  labs(title = "Rank 60~79") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) +
  coord_flip() -> p4


rank_range = c(80:99)
rank_pattern <- paste(rank_range, collapse = "|")
rank_pattern

plot.data %>% 
  filter(str_detect(Rank, rank_pattern)) %>% #head
  mutate(name = factor(name, levels = c("Missmatch", "Match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = value), 
            position = position_fill(vjust = 0.5)) + 
  labs(title = "Rank 80~99") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) +
  coord_flip() -> p5

p6 <- get_legend(plot.data %>% 
                   filter(str_detect(Rank, rank_pattern)) %>% #head
                   mutate(name = factor(name, levels = c("Missmatch", "Match"))) %>%  # Fill 순서 반대로 설정
                   ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
                   geom_bar(stat='identity',position = 'fill')+
                   guides(fill = guide_legend(reverse=T)) +
                   theme(legend.position = "bottom",
                         legend.title = element_blank()))



cowplot::plot_grid(p1,p2,p3,p4,p5,nrow = 1,label_x = "Percentage") -> p11
gridExtra::grid.arrange(p11,bottom = "Percentage",top="# of STRs",left = "Repeat Unit") -> p11

cowplot::plot_grid(p11,p6,nrow = 2,rel_heights = c(1,0.1))





head(str_RU_count)
head(str_RU_count_bymatch)

str_match_score_freqeuncy %>% mutate(match = ifelse(concordance == 1,"Match","Missmatch")) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>% #head()
  group_by(RU.length,match) %>%
  count(RU.length) %>% #head()
  pivot_wider(names_from = match,values_from = n) -> str_RU.length_count_bymatch

head(str_RU_count)
head(str_RU.length_count_bymatch)

str_RU.length_count_bymatch %>% mutate(n = Match + Missmatch) %>% 
  pivot_longer(cols = c("Match","Missmatch")) %>% #head()
  mutate(name = factor(name, levels = c("Missmatch", "Match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=RU.length,y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill')  + 
  geom_text(aes(label = value), 
            position = position_fill(vjust = 0.5)) + 
  scale_y_continuous(name = "Percentage",labels = percent) + 
    xlab("Repeat Unit Length") + 
  theme(legend.position = "bottom",legend.title = element_blank(),
      axis.text = element_text(size = 10)) + 
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip(clip = 'off')  + 
    scale_x_reverse(breaks = 2:24)

str_match_score_freqeuncy %>% #head()
  mutate(match = ifelse(concordance == 1,"ALL match",
                        ifelse(concordance > 0.5, "0.5<x<1",
                               ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  mutate(match = factor(match,levels=c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch"))) %>%
  mutate(RU.length = str_length(MOTIFS)) %>% #head()
  group_by(RU.length,match) %>%
  count(RU.length) %>% #head()
  pivot_wider(names_from = match,values_from = n) -> str_RU.length_count_bymatch

head(str_RU_count)
head(str_RU.length_count_bymatch)

str_RU.length_count_bymatch %>% #mutate(add_text = paste0(`0<x=<0.5`,"   ",`ALL missmatch`)) %>%
  mutate(add_text = paste0(`0<x=<0.5`,"        ")) %>%
  mutate(add_text = str_replace(add_text,"NA","")) %>% #head()
#  mutate(add_text = ifelse())
  pivot_longer(cols = c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch")) %>% #head()
  mutate(name = factor(name, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=RU.length,y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill')  + 
  geom_text(aes(label = ifelse(name %in% c("ALL match"),value,""), y = 0),
            hjust = 0, nudge_x = 0.05) + 
  geom_text(aes(label = ifelse(name %in% c("0.5<x<1"),value,"")),
            position = position_fill(vjust = 0)) + 
  geom_text(aes(label = ifelse(name %in% c("ALL missmatch"),value,""),fontface = "bold",y = 1), 
            hjust = 0, nudge_x = 0.05,color = "darkred",) + 
  geom_text(aes(label = ifelse(name %in% c("0<x=<0.5"),add_text,"")), 
            position = position_fill(vjust = 0),color='darkgreen') + 
  scale_y_continuous(name = "Percentage",labels = percent) + 
  xlab("Repeat Unit Length (bp)") + 
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.text = element_text(size = 10)) + 
  guides(fill = guide_legend(reverse=T)) + 
  #guides(fill = guide_legend(reverse=T,nrow=2,byrow=TRUE)) + 
  #guides(fill=guide_legend(nrow=2,byrow=TRUE))
  coord_flip(clip = 'off')  + 
  scale_x_reverse(breaks = 2:24)



  
  
######### concordacne 범위
str_match_score_freqeuncy <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt")
str_RU_count<- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.RUcount.txt",header =T)
str_qc_pass_common <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.txt",header = T)
head(str_RU_count)
colnames(str_RU_count)[1] <- "MOTIFS"

#head(str_match_score_freqeuncy)
str_match_score_freqeuncy %>% #mutate(match = ifelse(concordance == 1,"Match","Missmatch")) %>% 
  mutate(match = ifelse(concordance == 1,"ALL match",
                                    ifelse(concordance > 0.5, "0.5<x<1",
                                           ifelse(concordance > 0,"0<x=<0.5","ALL missmatch")))) %>%
  mutate(match = factor(match,levels=c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch"))) %>%
  group_by(MOTIFS,match) %>%
  count(MOTIFS) %>% #head()
  pivot_wider(names_from = match,values_from = n) -> str_RU_count_bymatch


head(str_RU_count)
head(str_RU_count_bymatch)

head(str_RU_count)

rank_min = 1
rank_max = 99


str_RU_count %>% left_join(str_RU_count_bymatch) %>% #head()
  select(-n) %>% #head()
  mutate(add_text = paste0(`0<x=<0.5`,"        ")) %>%
  mutate(add_text = str_replace(add_text,"NA","")) %>% #head()
  pivot_longer(cols = c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch")) %>% #head()
  mutate(Rank = ifelse(Rank %in% c(rank_min:rank_max),Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Rank %in% c(rank_min:rank_max),MOTIFS,"Other")) %>% na.omit() %>% #head()
  group_by(MOTIFS,Rank,name) %>% 
  #summarise(value = sum(value)) %>% 
  left_join(str_RU_count %>% select(-Rank)) %>% #head()\
  filter(MOTIFS != "Other")-> plot.data

head(plot.data)

str_RU_count %>% left_join(str_RU_count_bymatch) %>% #head()
  pivot_longer(cols = c("ALL match","0.5<x<1",0,"0<x=<0.5","ALL missmatch")) %>% #head()
  mutate(Rank = ifelse(Rank %in% c(rank_min:rank_max),Rank,"Other")) %>%
  mutate(MOTIFS = ifelse(Rank %in% c(rank_min:rank_max),MOTIFS,"Other")) %>% na.omit() %>% 
  filter(Rank == "Other") %>%
  group_by(MOTIFS,Rank,name) %>%
  summarise(value = sum(value)) %>% ungroup() %>% 
  mutate(add_text = paste0(value,"        ")) %>%
  mutate(n = sum(value)) -> plot.data.other

head(plot.data)
head(plot.data.other)

#plot.data %>% dim()
as.character(c(1:19))

#head(plot.data %>% filter(Rank %in% c(1:20-1)))
head(plot.data)
plot.data %>% filter(Rank %in% as.character(c(1:19))) %>% rbind(plot.data.other)%>% #filter(MOTIFS == "Other")
  mutate(name = factor(name, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  #geom_text(aes(label = ifelse(name %in% c("0.5<x<1","ALL match"),value,"")),
  geom_text(aes(label = ifelse(name %in% c("ALL match"),value,""), y = 0),
            hjust = 0) + 
  geom_text(aes(label = ifelse(name %in% c("0.5<x<1"),value,""), y = 0.5),
            hjust = 0) + 
  
  geom_text(aes(label = ifelse(name %in% c("ALL missmatch"),value,""),fontface = "bold",y = 1), 
            hjust = 0,color = "darkred",) + 
  geom_text(aes(label = ifelse(name %in% c("0<x=<0.5"),add_text,"")), 
            position = position_fill(vjust = 0),color='darkgreen') + 
  scale_y_continuous(labels = percent) + 
    labs(title = "Other, Rank 1~19") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p1
  
p1

rank_range = as.character(c(20:39))
plot.data %>% filter(Rank %in% rank_range) %>% #rbind(plot.data.other)%>% #filter(MOTIFS == "Other")
    mutate(name = factor(name, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = ifelse(name %in% c("ALL match"),value,""), y = 0),
            hjust = 0, nudge_x = 0.05) + 
  geom_text(aes(label = ifelse(name %in% c("0.5<x<1"),value,"")),
            position = position_fill(vjust = 0)) + 
  geom_text(aes(label = ifelse(name %in% c("ALL missmatch"),value,""),fontface = "bold",y = 1), 
            hjust = 0, nudge_x = 0.05,color = "darkred",) + 
  geom_text(aes(label = ifelse(name %in% c("0<x=<0.5"),add_text,"")), 
            position = position_fill(vjust = 0),color='darkgreen') + 
  scale_y_continuous(labels = percent) + 
  labs(title = "Other, Rank 20~39") +
  theme(legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p2
p2

p1
p2
#max(plot.data[plot.data$Rank %in% rank_range,]$value)
#plot.data %>% filter(MOTIFS == "AAAAC")

rank_range = as.character(c(40:59))
rank_range
plot.data %>% filter(Rank %in% rank_range) %>% #rbind(plot.data.other)%>% #filter(MOTIFS == "Other")
  mutate(name = factor(name, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = ifelse(name %in% c("ALL match"),value,""), y = 0),
            hjust = 0, nudge_x = 0.05) + 
  geom_text(aes(label = ifelse(name %in% c("0.5<x<1"),value,"")),
            position = position_fill(vjust = 0)) + 
  geom_text(aes(label = ifelse(name %in% c("ALL missmatch"),value,""),fontface = "bold",y = 1), 
            hjust = 0, nudge_x = 0.05,color = "darkred",) + 
  geom_text(aes(label = ifelse(name %in% c("0<x=<0.5"),add_text,"")), 
            position = position_fill(vjust = 0),color='darkgreen') + 
  scale_y_continuous(labels = percent) + 
  labs(title = "Other, Rank 40~59") +
  theme(legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p3
p3
p2

rank_range = as.character(c(60:79))
plot.data %>% filter(Rank %in% rank_range) %>% #rbind(plot.data.other)%>% #filter(MOTIFS == "Other")
  mutate(name = factor(name, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = ifelse(name %in% c("ALL match"),value,""), y = 0),
            hjust = 0, nudge_x = 0.05) + 
  geom_text(aes(label = ifelse(name %in% c("0.5<x<1"),value,"")),
            position = position_fill(vjust = 0)) + 
  geom_text(aes(label = ifelse(name %in% c("ALL missmatch"),value,""),fontface = "bold",y = 1), 
            hjust = 0, nudge_x = 0.05,color = "darkred",) + 
  geom_text(aes(label = ifelse(name %in% c("0<x=<0.5"),add_text,"")), 
            position = position_fill(vjust = 0),color='darkgreen') + 
  scale_y_continuous(labels = percent) + 
  labs(title = "Other, Rank 60~79") +
  theme(legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p4
p4

rank_range = as.character(c(80:99))
plot.data %>% filter(Rank %in% rank_range) %>% #rbind(plot.data.other)%>% #filter(MOTIFS == "Other")
  mutate(name = factor(name, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  geom_text(aes(label = ifelse(name %in% c("ALL match"),value,""), y = 0),
            hjust = 0, nudge_x = 0.05) + 
  geom_text(aes(label = ifelse(name %in% c("0.5<x<1"),value,"")),
            position = position_fill(vjust = 0)) + 
  geom_text(aes(label = ifelse(name %in% c("ALL missmatch"),value,""),fontface = "bold",y = 1), 
            hjust = 0, nudge_x = 0.05,color = "darkred",) + 
  geom_text(aes(label = ifelse(name %in% c("0<x=<0.5"),add_text,"")), 
            position = position_fill(vjust = 0),color='darkgreen') + 
  scale_y_continuous(labels = percent) + 
  labs(title = "Other, Rank 80~99") +
  theme(legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p5
p5

plot.data %>% filter(Rank %in% rank_range) %>% #rbind(plot.data.other)%>% #filter(MOTIFS == "Other")
  mutate(name = factor(name, levels = c("ALL missmatch","0<x=<0.5","0.5<x<1","ALL match"))) %>%  # Fill 순서 반대로 설정
  ggplot(aes(x=reorder(MOTIFS,n),y=value,fill=name)) + 
  geom_bar(stat='identity',position = 'fill') + 
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(fill = guide_legend(reverse=T)) + 
  coord_flip() -> p6


p6 <- get_legend(p6)
p6

cowplot::plot_grid(p1,p2,p3,p4,p5,nrow = 1,label_x = "Percentage") -> p11
gridExtra::grid.arrange(p11,bottom = "Percentage",top="# of STRs",left = "Repeat Unit") -> p11

cowplot::plot_grid(p11,p6,nrow = 2,rel_heights = c(1,0.1))



##### RU type by non match
str_match_score_freqeuncy <- read_table("~/Desktop/KU/@research/STR/02.compare/STR_type/str_match_score_frequency.txt")
str_RU_count<- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.RUcount.txt",header =T)
str_qc_pass_common <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.qcpass.common.txt",header = T)
head(str_RU_count)
colnames(str_RU_count)[1] <- "MOTIFS"

head(str_RU_count)
head(str_qc_pass_common)
head(str_match_score_freqeuncy)

str_match_score_freqeuncy %>% filter(concordance == 0) %>% 
  count(MOTIFS) %>% filter(n != 1) %>%
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(pct=n/sum(n)*100) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(MOTIFS = paste0(MOTIFS,'(',RU.length,')') ) %>%
  arrange(n, -RU.length) %>% #filter(n==2)
  mutate(MOTIFS = factor(MOTIFS, levels = unique(MOTIFS))) %>%
  ggplot(aes(x=MOTIFS,y=n,fill=factor(RU.length))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(n)),
           position=position_stack(vjust=0.5)) + 
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "None",
        plot.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
  axis.text = element_text(size = 10)) + 
  guides(fill = guide_legend(reverse=F)) + 
  coord_flip(clip = 'off') -> p1

p1  

str_match_score_freqeuncy %>% filter(concordance == 0) %>% 
  count(MOTIFS) %>% filter(n == 1) %>%
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(pct=n/sum(n)*100) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(MOTIFS = paste0(MOTIFS,'(',RU.length,')') ) %>%
  arrange(n, -RU.length) %>% #filter(n==2)
  mutate(MOTIFS = factor(MOTIFS, levels = unique(MOTIFS))) %>%
  ggplot(aes(x=MOTIFS,y=n,fill=factor(RU.length))) + 
  geom_bar(stat='identity') + 
  #geom_text(aes(label=paste0(n)),
            #position=position_stack(vjust=0.5)) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "None",
        axis.text = element_text(size = 10)) + 
  guides(fill = guide_legend(reverse=F)) + 
  scale_y_continuous(breaks = c(0, 1)) + 
  coord_flip(clip = 'off') -> p2
p2

str_match_score_freqeuncy %>% filter(concordance == 0) %>% 
  count(MOTIFS) %>% arrange(-n) %>% #filter(n == 1) %>%
  group_by(MOTIFS) %>% #head()
  summarise(n = sum(n)) %>% #dim()
  mutate(RU.length = str_length(MOTIFS)) %>%
  ggplot(aes(x=reorder(MOTIFS,n),y=n,fill=factor(RU.length))) + 
  geom_bar(stat='identity') + 
  #labs(caption = "111") + 
  #guide_legend(title = "hit") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 10)) + 
  guides(fill = guide_legend(reverse=F,title="RU length")) -> p3

p3 <- get_legend(p3)



cowplot::plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(5,3,1)) -> p11

gridExtra::grid.arrange(p11, bottom = "# of STRs",left = "Repeat Unit (RU length)")


