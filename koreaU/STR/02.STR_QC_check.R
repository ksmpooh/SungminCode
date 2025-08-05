library(tidyverse)
library(data.table)
library(ggbreak)
library(ggpubr)
library(ggplot2)
library(cowplot)

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/Quality_check/")

## EH quality check

flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  tmp$ID <- str_replace(i,".EH.qcmt.txt","")
  df <- rbind(df,tmp)
}
#df <- rbind(df,tmp)
head(df)
colnames(df) <- c("chr","pos","STR_ID","RU","ALT","FILTER","GT","type","ADSP","ADFL","ADIR","ID")

df %>% select(ID,STR_ID,RU,FILTER,type) -> df_qc
head(df_qc)

df_qc %>% group_by(ID) %>% count(FILTER) -> df_filter
head(df_filter)
df_filter %>% ggplot(aes(x=ID,y=n,fill=FILTER)) +
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of STR") + 
  guides(fill=guide_legend(title="QV")) + 
  theme(#legend.position = "bottom",
        axis.text.x = element_blank())

head(df_filter)
df_filter %>% filter(FILTER == "LowDepth") %>%
  ggplot(aes(x=ID,y=n)) +
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of LowDepth STR") + 
  theme(legend.position = "bottom",
        axis.text.x = element_blank())


head(df_qc)
df_qc %>% select(ID,FILTER,type) %>% group_by(ID,type) %>% count(FILTER) -> df_filter_withtype

head(df_filter_withtype)
df_filter_withtype %>% mutate(type = ifelse(type == "./.","Missing value",type)) %>%
  ggplot(aes(x=ID,y=n,fill=type)) + 
  geom_bar(position="stack", stat="identity")

table(df_filter_withtype$type)

df_filter_withtype %>% mutate(type = ifelse(type == "./.","Missing value",type)) %>% filter(FILTER == "LowDepth") %>%
  ggplot(aes(x=ID,y=n,fill=type)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of LowDepth STR") + 
  guides(fill=guide_legend(title="Read type")) + 
  #  scale_fill_manual(values = c("Missing value" = "grey")) + 
  theme(axis.text.x = element_blank())

library(unikn)

palette(usecol(pal_unikn, n = 9))

palette(rainbow(n = 9)) -> a
c("grey",rev(a)) -> a
a

unique(df_filter_withtype$type)[2:10] %>% sort() -> b
c("Missing",b) -> b
b
a
df_filter_withtype %>% mutate(type = ifelse(type == "./.","Missing",type)) %>% filter(FILTER == "LowDepth") %>%
  ggplot(aes(x=ID,y=n,fill=factor(type, levels = b))) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of LowDepth STR") + 
  guides(fill=guide_legend(title="Read type")) + 
  #scale_fill_manual(values = a,breaks = b) +
  scale_fill_manual(values = a) + 
  theme(axis.text.x = element_blank())

head(df_filter_withtype)
df_filter_withtype %>% mutate(type = ifelse(type == "./.","Missing",type)) %>% filter(FILTER != "LowDepth") -> df_filter_withtype_pass
df_filter_withtype_pass

df_filter_withtype %>% mutate(type = ifelse(type == "./.","Missing",type)) %>% filter(FILTER != "LowDepth") %>%
  ggplot(aes(x=ID,y=n,fill=factor(type, levels = b))) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of STR") + 
  guides(fill=guide_legend(title="Read type")) + 
  #scale_fill_manual(values = a,breaks = b) +
  scale_fill_manual(values = a[2:10]) + 
  theme(axis.text.x = element_blank()) + 
  coord_cartesian(ylim = c(317000,321000))





#### TGRT
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/Quality_check/")

## TGRT quality check
#system.time(tmp <- read_table(i))
#system.time(tmp <- fread(i))

head(tmp)
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
count = 0
for (i in flist) {
  #tmp <- read_table(i)
  count = count + 1
  print(count)
  tmp <- fread(i)
  tmp$ID <- str_replace(str_replace(i,".pbmm2_hg38_withunmapped_trgt_genotype_gangstr.sorted.MC_AP_AM.txt",""),"_sorted_trgt_genotype_gangstr.sorted.MC_AP_AM.txt","")
  df <- rbind(df,tmp)
}
#df <- rbind(df,tmp)
head(df)
#colnames(df) <- c("chr","pos","STR_ID","RU","ALT","FILTER","GT","type","ADSP","ADFL","ADIR","ID")

df %>% select(-CHROM,-POS,-STRUC,-ALT) -> df_mc_ap_am
head(df_mc_ap_am)

df_mc_ap_am %>% filter(!(AP %in% c(".","NA"))) %>% filter(MC != 0) %>% mutate(AP = as.numeric(AP)) -> df_mc_ap_am
#df_mc_ap_am %>% filter(!(AP %in% c(".","NA"))) %>% filter(MC == 0) -> trgt_check
df_mc_ap_am %>% filter(str_detect(MC,"_")) -> trgt_check
#head(trgt_check)
#write.table(trgt_check,"~/Desktop/KU/@research/STR/figure/trgt_complex_str_basic.info.txt",col.names = T,row.names = F,quote = F,sep = "\t")
head(df_mc_ap_am)

df_mc_ap_am %>% ggplot(aes(x=ID,y=AP)) + 
  geom_violin() + 
  xlab("Sample") +
  ylab("AP") + 
  coord_cartesian(ylim = c(0.95,1)) + 
  theme(axis.text.x = element_blank())


df_mc_ap_am %>% filter(AP == 1) %>% group_by(ID) %>% count(ID) %>% #summary(n)
  ggplot(aes(x=ID,y=n)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Sample") +
  ylab("# of (AP = 1)") + 
  coord_cartesian(ylim = c(563000,566500)) + 
  theme(axis.text.x = element_blank())



df_mc_ap_am %>% filter(AP == 1) %>% group_by(TRID) %>% count(TRID) -> a
df_mc_ap_am %>% filter(AP >= 0.8) %>% group_by(TRID) %>% count(TRID) -> a

head(a)
a %>% filter(n == 132) %>% dim()


df %>% select(-CHROM,-POS,-STRUC,-ALT)  %>% filter(!(AP %in% c(".","NA"))) %>% mutate(AP = as.numeric(AP)) -> df_mc_ap_am

head(df_mc_ap_am)
#df_mc_ap_am %>% pivot_wider(names_from = Allele,values_from = MC) -> df_mc_ap_am_re

df_mc_ap_am %>% group_by(ID,TRID) %>% summarise(allele_count = n(),mean_AP = mean(AP)) -> df_mc_ap_am_check

head(df_mc_ap_am_check)
df_mc_ap_am_check %>% filter(allele_count == 2,mean_AP== 1) %>% dim()
#df_mc_ap_am_check %>% filter(mean_AP == 1) %>% dim()

head(df)
#table(df$TRID) %>% dim


#### distribution
load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.info.RData")
load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.info.RData")

trgt <- tgrt
tgrt <- NULL
head(tgrt)
head(eh)
head(trgt)

eh %>% select()
eh %>% filter(FILTER == "PASS") %>% select(-type,-ADSP,-ADFL,-ADIR)

eh %>% select(-type,-ADSP,-ADFL,-ADIR) %>% filter(FILTER == "PASS") %>%
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>%
  mutate(ALT = ifelse(ALT == ".",0,ALT)) %>% mutate(GT = str_split_fixed(GT,"/",2),ALT= str_split_fixed(ALT,",",2)) -> eh_pass_gt

eh_pass_gt %>% mutate(STR1 = ifelse(GT[,1] == "1",ALT[,1],ifelse(GT[,1] == "2",ALT[,2],REF))) %>%
  mutate(STR2 = ifelse(GT[,2] == "1",ALT[,1],ifelse(GT[,2] == "2",ALT[,2],REF))) %>% select(-ALT,-GT,-REF) %>% mutate(STR1 = as.integer(STR1),STR2 = as.integer(STR2)) -> eh_pass_gt


eh %>% filter(FILTER == "PASS") %>% select(ID,STR_ID,type) %>% mutate(type1 = str_split_fixed(type,"/",2)[,1],type2 = str_split_fixed(type,"/",2)[,2]) %>% 
  select(-type) -> eh_pass_type

head(eh_pass_type)
head(eh_pass_gt)

eh_pass_gt %>% select(-FILTER) %>% left_join(eh_pass_type) -> eh_pass

head(eh_pass)
eh_pass %>% mutate(EU_STR1 = ifelse(STR1 < STR2, STR1,STR2),EU_STR2 = ifelse(STR1 < STR2, STR2,STR1)) %>%
  mutate(EU_type1 = ifelse(STR1 < STR2, type1,type2),EU_type2 = ifelse(STR1 < STR2, type2,type1)) %>% 
  select(-STR1,-STR2,-type1,-type2)-> eh_pass

head(eh_pass)



head(trgt)
trgt %>% filter(AP >=0.8) %>% dim()
trgt %>% filter(MC == 0) %>% head(10)

trgt %>% filter(!(AP %in% c(".","NA"))) %>% filter(MC != 0) %>% select(-AM) %>% #dim()
  mutate(AP = as.numeric(AP)) %>% filter(AP >=0.8) -> trgt_pass



head(trgt_pass)
head(eh_pass)


common_str <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/TRID_common_EHpass_TRGTupper0.8.txt",header = T)


head(ref)
head(common_str)

head(trgt_pass)
head(eh_pass)


trgt_pass %>% filter(TRID %in% common_str$TRID) -> trgt_pass_common
eh_pass %>% filter(STR_ID %in% common_str$TRID) -> eh_pass_common

trgt_pass_common %>% mutate(ID = str_replace_all(ID,".merge","")) -> trgt_pass_common

#colnames(trgt_pass_common) <- c("")
#save(eh_pass_common,file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.PASS.intersect.RData")
#save(trgt_pass_common,file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.upper0.8.intersect.RData")

#####
load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.PASS.intersect.RData")
load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.upper0.8.intersect.RData")

ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)
common_str <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/TRID_common_EHpass_TRGTupper0.8.txt",header = T)

head(common_str)
head(trgt_pass_common)
head(eh_pass_common)
#trgt_pass_common %>% filter(ID == "NIH23F1013274") %>% select(TRID) %>% unique() %>% dim()
trgt_pass_common %>% filter(ID == "NIH23F1013274") %>%
  select(TRID,MOTIFS) %>% count(MOTIFS)

head(eh_pass_common)

eh_pass_common %>% filter(ID == "NIH20N2000078") %>%
  select(STR_ID,RU) %>% count(RU) %>% arrange(-n) %>% mutate(Rank = rank(-n)) %>% 
  mutate(STR_MOTIFs = ifelse(Rank == c(1:10),RU,"other")) %>% #head()
  group_by(STR_MOTIFs) %>%
  summarise(n = sum(n)) %>%
  mutate(pct=n/sum(n)*100) %>% #head()
  ggplot(aes(x=reorder(STR_MOTIFs,-n),y=n,fill=reorder(STR_MOTIFs,-n))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(round(pct,1), '%')),
              position=position_stack(vjust=0.99)) + 
#  scale_y_cut(breaks=c(50000, 160000), scales=c(0)) +
  scale_y_break(c(50000, 150000)) +
  ylab("Count") + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 10))

eh_pass_common %>% filter(ID == "NIH20N2000078") %>%
  select(STR_ID,RU) %>% count(RU) %>% arrange(-n) %>% mutate(Rank = rank(-n)) -> str_RU_count
#  filter(Rank %in% c(1:10)) -> str_RU_count
#  mutate(STR_MOTIFs = ifelse(Rank == c(1:10),RU,"other")) %>% #head()
str_RU_count %>% filter(n == 1) %>% dim()

trans <- function(x){pmin(x,500) + 0.05*pmax(x-500,0)}
yticks <- c(0,100,200,400,600,800,4900,5000)

str_RU_count %>% mutate(new_n = ifelse(n %in% c(1:10),n,">10")) %>% 
  #group_by(RU) %>%
  count(new_n) %>% #head()
  ggplot(aes(x=factor(new_n,c("1","2","3","4","5","6","7","8","9","10",">10")),y=n,fill=new_n)) + 
  geom_bar(stat='identity') +
  scale_y_break(c(700, 4900),ticklabels = c(1,2,3,4950)) + 
#  ylim(c(0,5000)) + 
  xlab("RU type count") + 
  ylab("Count") + 
  theme(#axis.title.x = element_blank(),
        legend.position = 'None',
        axis.text = element_text(size = 10))


####
head(str_match_score_freqeuncy)
head(ref)
common_str <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/TRID_common_EHpass_TRGTupper0.8.txt",header = T)
ru <- read_table("~/Desktop/KU/@research/STR/eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed",col_names = F)
head(ru)
head(common_str)
head(ru)
ru %>% mutate(MOTIFS = str_split_fixed(X4,";",3)[,2]) %>% 
  mutate(MOTIFS = str_split_fixed(MOTIFS,"=",2)[,2]) %>% #head()
  mutate(ID = str_split_fixed(X4,";",3)[,1]) %>% 
  mutate(ID = str_split_fixed(ID,"=",2)[,2]) %>% #head()
  select(X1,X2,X3,MOTIFS,ID) -> ru

head(ru)
head(common_str)
colnames(ru) <- c("chrom","start","end","MOTIFS","STR")
ru %>% mutate(STR_region = paste0(chrom,"_",start,"_",end)) %>% filter(STR %in% common_str$TRID)-> ru

#write.table(ru,"/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.IDregion.txt",row.names = F,col.names = T,quote = F,sep = "\t")
####
####

#ru,"/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.IDregion.txt"
ru <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.IDregion.txt")
head(ru)
ru.scale = c(seq(1,14),"15+")
ru.scale
ru %>% select(MOTIFS) %>% mutate(RU.length = str_length(MOTIFS)) %>%
  count(RU.length) %>% mutate(RU.length = ifelse(RU.length >= 15,"15+",RU.length)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  labs(x="RU length",y="# of STRs") -> p1

ru %>% select(MOTIFS) %>% mutate(RU.length = str_length(MOTIFS)) %>%
  count(RU.length) %>% mutate(RU.length = ifelse(RU.length >= 15,"15+",RU.length)) %>%
  ggplot(aes(x=factor(RU.length,levels=ru.scale),y=n)) + 
  geom_bar(stat = 'identity') + 
  xlim(c(seq(7,14),"15+")) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())-> p2


combined_plot <- ggdraw() +
  draw_plot(p1, 0, 0, 1, 1) +           # p1�� ���� ũ��� ��ġ
  draw_plot(p2, 0.5, 0.5, 0.5, 0.5)     # p2�� ������ ���� �۰� ��ġ

# ��� ���
print(combined_plot)

head(ru)

ru %>%
  select(STR,MOTIFS) %>% count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) %>% 
  mutate(STR_MOTIFs = ifelse(Rank %in% c(1:10),MOTIFS,"other")) %>% #head()
  group_by(STR_MOTIFs) %>%
  summarise(n = sum(n)) %>%
  mutate(pct=n/sum(n)*100) %>% #head()
  ggplot(aes(x=reorder(STR_MOTIFs,-n),y=n,fill=reorder(STR_MOTIFs,-n))) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=paste0(round(pct,1), '%')),
            position=position_stack(vjust=0.99)) + 
#    scale_y_cut(breaks=c(50000, 160000), scales=c(0)) +
  scale_y_break(c(50000, 150000)) +
  labs(y="# of STRs") +
  theme(axis.title.x = element_blank(),
        #legend.title = element_blank(),
        legend.position = "None",
        axis.text = element_text(size = 10)) -> p1

p1

head(ru)
ru %>%
  select(STR,MOTIFS) %>% count(MOTIFS) %>% arrange(-n) %>% mutate(Rank = rank(-n)) -> str_RU_count

trans <- function(x){pmin(x,500) + 0.05*pmax(x-500,0)}
yticks <- c(0,100,200,400,600,800,4900,5000)

str_RU_count %>% mutate(new_n = ifelse(n %in% c(1:10),n,">10")) %>% 
  #group_by(RU) %>%
  count(new_n) %>% #head()
  ggplot(aes(x=factor(new_n,c("1","2","3","4","5","6","7","8","9","10",">10")),y=n,fill=new_n)) + 
  geom_bar(stat='identity') +
    scale_y_break(c(1000, 4700),ticklabels = c(1,2,3,4900,5000)) + 

  labs(x="# of RU types",y="# of STRs")  + 
  theme(#axis.title.x = element_blank(),
    legend.position = 'None',
    #axis.y.
    axis.text = element_text(size = 10)) -> p2
p2
#p2

combined_plot <- ggdraw() +
  draw_plot(p1, 0, 0, 1, 1) +           # p1�� ���� ũ��� ��ġ
  draw_plot(p2, 0.5, 0.5, 0.5, 0.5)     # p2�� ������ ���� �۰� ��ġ

# ��� ���
print(combined_plot)


#install.packages("tidyverse")
library(tidyverse)
##### STR type
ref_db <- readxl::read_xlsx("~/Desktop/KCDC/paper/STR/2023_guo_sup_2.xlsx",sheet = 1)
head(ref_db)
ref_db %>% select(chr,`STR start position (GRCh38)`,`STR end position (GRCh38)`,motif,`Reference tract length`,`Genomic annotation`,`Distance to nearest TSS`) -> ref_db

colnames(ref_db) <- c("chrom","start","end","MOTIFS","Referencetractlength","anotation","DistancetoTSS")
ref_db %>% mutate(STR_region = paste0(chrom,"_",start,"_",end)) -> ref_db

ru <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.IDregion.txt")
dim(ru)

head(ru)
head(ref_db)

library(plotly)
library(ggrepel)

ru %>% left_join(ref_db) %>% mutate(anotation = ifelse(is.na(anotation),"chrX",anotation)) %>% #head()
  group_by(anotation) %>% #head()
  count() %>% ungroup() %>%
  mutate(perc = n/sum(n)) %>%
  arrange(perc) %>% #-> a
  mutate(labels = scales::percent(perc)) %>% #head
  ggplot(aes(1, perc, fill = anotation)) +
  geom_col(color = 'black', 
           position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  geom_text_repel(aes(x = 1.4, y = , label = label), 
                  nudge_x = .3, 
                  segment.size = .7, 
                  show.legend = FALSE) +
  coord_polar('y') +
  theme_void()

  






dup <- read_table('~/Desktop/KU/@research/STR/db/genomicSuperDups.txt.gz',col_names = F)
head(dup)
dup[1:5,]

###


                                       

str_RU_count %>% filter(Rank %in% c(1:10)) -> str_RU_count_10rank
str_RU_count %>% filter(n >= 10) %>% dim()

head(trgt_pass_common)
head(eh_pass_common)

trgt_pass_common %>% select(ID,TRID,MOTIFS,Allele,MC) %>% #head()
  pivot_wider(names_from = Allele,values_from = MC) %>% 
  mutate(New_MC = (allele_1+allele_2)/2) %>% group_by(MOTIFS) %>%
  summarise(mean_MC_TRGT = mean(New_MC)) -> trgt_pass_common_meanMC
  

eh_pass_common %>% select(ID,STR_ID,RU,EU_STR1,EU_STR2) %>%
  mutate(New_MC = (EU_STR1+EU_STR2)/2) %>% group_by(RU) %>% 
  summarise(mean_MC_EH = mean(New_MC)) -> eh_pass_common_meanMC

head(trgt_pass_common_meanMC)
head(eh_pass_common_meanMC)

trgt_pass_common_meanMC %>% rename(RU = MOTIFS) %>% left_join(eh_pass_common_meanMC) %>% #head()
  left_join(str_RU_count_10rank) %>% mutate(STR = ifelse(is.na(Rank),NA,RU)) %>% #head()
  mutate(diff = mean_MC_TRGT - mean_MC_EH) %>% #ggplot(aes(y=diff,x=RU)) + geom_point() 
  mutate(STR_mean_diff = ifelse(abs(diff) >1,">1","1>=")) %>%
  ggpubr::ggscatter(.,x='mean_MC_TRGT',y="mean_MC_EH",fill = "STR_mean_diff",color = "STR_mean_diff",#,color='Rank',
          add = "reg.line",
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          #xlim = c(0,0.5),
          #ylim = c(0,0.5),
          #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = "TRGT",
          ylab = "EH") + 
  geom_abline(slope = 1,linetype="dashed") + 
  theme(legend.position = "right")


trgt_pass_common_meanMC %>% rename(RU = MOTIFS) %>% left_join(eh_pass_common_meanMC) %>% #head()
  mutate(diff = mean_MC_TRGT - mean_MC_EH) %>% #ggplot(aes(y=diff,x=RU)) + geom_point() 
  mutate(RU.length = str_length(RU)) %>%
  mutate(STR_mean_diff = ifelse(abs(diff) >1,">1","1>=")) %>%
  #ggpubr::ggscatter(.,x='mean_MC_TRGT',y="mean_MC_EH",fill = "STR_mean_diff",color = "STR_mean_diff",size = "RU.length",
  ggscatter(.,x='mean_MC_TRGT',y="mean_MC_EH",color = "RU.length",#,shape = "STR_mean_diff",
                    add = "reg.line",
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #xlim = c(0,0.5),
                    #ylim = c(0,0.5),
                    #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
                    cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                    xlab = "TRGT",
                    ylab = "EH") + 
  geom_abline(slope = 1,linetype="dashed") + 
  gradient_color(c("white","red")) + 
  theme(legend.position = "right")
  

trgt_pass_common_meanMC %>% rename(RU = MOTIFS) %>% left_join(eh_pass_common_meanMC) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",RU))/str_length(RU)) %>%
  mutate(RU.length = str_length(RU)) %>%
  #ggscatter(.,x='mean_MC_TRGT',y="mean_MC_EH",color = "GC",size='RU.length',
  ggscatter(.,x='mean_MC_TRGT',y="mean_MC_EH",size = "GC",color='RU.length',
            add = "reg.line",
            conf.int = TRUE, # Add confidence intervanamel
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            #xlim = c(0,0.5),
            #ylim = c(0,0.5),
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n"),
            xlab = "TRGT",
            ylab = "EH") + 
  geom_abline(slope = 1,linetype="dashed") + 
  gradient_color(c("white","blue")) + 
#  guides(color=guide_legend(title="GC contents of RU")) + 
  theme(legend.position = "right")


head(eh_pass_common)
eh_pass_common %>% count(ID) %>% count(n)
head(trgt_pass_common)

trgt_pass_common %>% count(ID) %>% count(n)

trgt_pass_common %>% select(ID,TRID,MOTIFS,Allele,MC) %>% #head()
  pivot_wider(names_from = Allele,values_from = MC) -> trgt_pass_common_gt


eh_pass_common %>% select(ID,STR_ID,RU,EU_STR1,EU_STR2) -> eh_pass_common_gt

head(eh_pass_common_gt)
head(trgt_pass_common_gt)


colnames(eh_pass_common_gt) <- c("EH_ID","STR","MOTIFS","EH_STR1","EH_STR2")
colnames(trgt_pass_common_gt) <- c("ID","STR","MOTIFS","TRGT_STR1","TRGT_STR2")

#trgt_pass_common_gt %>% filter(is.na(TRGT_STR2))
ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)

ref %>% select(Revio,Illumina) -> ref
head(ref)
colnames(ref) <- c("ID","EH_ID")
table(eh_pass_common_gt$EH_ID %in% trgt_pass_common_gt$ID)


eh_pass_common_gt %>% left_join(ref) %>%  select(EH_ID,ID) %>% unique()-> id1
trgt_pass_common_gt %>% select(ID) %>% unique() -> id2

head(id1)
head(id2)

id1 %>% filter(!(ID %in% id2$ID))
id2 %>% filter(!(ID %in% id1$ID))


eh_pass_common_gt %>% left_join(ref) %>% #filter(ID %in% trgt_pass_common_gt$ID) %>% count(ID)
  select(-EH_ID) %>% left_join(trgt_pass_common_gt) -> df

head(df)
df %>% filter(is.na(TRGT_STR1))
df %>% select(1:7) %>% rename(tr_st1 = TRGT_STR1,tr_st2 = TRGT_STR2) %>% #head()
  mutate(TRGT_STR1 = ifelse(tr_st1 < tr_st2,tr_st1,tr_st2),TRGT_STR2 = ifelse(tr_st1 < tr_st2,tr_st2,tr_st1)) %>% select(-tr_st1,-tr_st2)-> df

head(df)
#df %>% count(EH_STR1 <= EH_STR2)
#df %>% count(TRGT_STR1 <= TRGT_STR2)
df %>% #mutate()
  mutate(STR1=TRGT_STR1-EH_STR1,STR2=TRGT_STR2-EH_STR2) %>% #head()
  #select(-EH_STR1,-EH_STR2,-TRGT_STR1,-TRGT_STR2) %>% 
  mutate(nRU_diff_mean = (STR1+STR2)/2) %>% select(-tr_st1,-tr_st2) -> df

head(df)

head(df)
df %>% count(EH_STR1 > EH_STR2)
df %>% count(TRGT_STR1 > TRGT_STR2)


######################### comapre
#write.table(df,"~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReaptnumber.with_dfiff.txt",col.names = T,row.names = F,quote = F,sep = "\t")
df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReaptnumber.with_dfiff.txt")
head(df)
df %>% count(ID) %>% count(n)
summary(df$nRU_diff_mean)

df %>% ggplot(aes(x=ID,y=nRU_diff_mean,fill=nRU_diff_mean)) + 
  geom_violin()+
  scale_y_break(c(50,300))+
  #scale_y_break(c(700, 4900),ticklabels = c(1,2,3,4950))
  theme(axis.text.x = element_blank())


df %>% filter(nRU_diff_mean != 0) %>%
  ggplot(aes(x=ID,y=nRU_diff_mean,fill=nRU_diff_mean)) + 
  geom_violin()+
  ylim(c(-1,1)) + 
  theme(axis.text.x = element_blank())



head(df)
head(df)

df %>% filter(STR1 != 0 | STR2 != 0) %>% 
  count(STR,MOTIFS)  -> b


df %>% filter(STR1 != 0,STR2 != 0) %>% 
  count(STR,MOTIFS)  -> c

df %>% filter(STR1 == 0,STR2 == 0)  %>% select(ID,STR,MOTIFS) %>% count(STR,MOTIFS) -> a

head(c)
head(a)
head(b)


#write.table(a,"~/Desktop/KU/@research/STR/02.compare/STR.longvsshort.match.count.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#a <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.longvsshort.match.count.txt")

a %>% count(n) %>% #head()
  ggplot(aes(x=n,y=nn)) + 
  geom_bar(stat='identity') + 
  scale_y_break(c(51000,140000),ticklabels = c(0,10000,20000,30000,40000,50000,145000))+
  xlab("# of Sample") +
  ylab("# of STR")


## str1 | str2 ���� �ϳ��� 0�ΰ�
b %>% count(n) %>% #head()
  ggplot(aes(x=n,y=nn)) + 
  geom_bar(stat='identity') + 
  scale_y_break(c(20000,45000),ticklabels = c(0,5000,10000,15000,20000,45000,50000))+
  xlab("# of Sample") +
  ylab("# of STR")

## str1, str2 �Ѵ� 0�ΰ�
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
c %>% filter(n==66) %>% dim()

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
        axis.text = element_text(size = 10))

  

head(a)
head(b)
head(str_RU_count)
head(all_sample_same_RUtype)

new_df %>% filter(TRGT_STR1 == 2)
#head(new_df)
head(df)
df %>% filter(STR %in% all_sample_same_RUtype$STR) %>% #head()
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(mean_repeat_count = (TRGT_STR1 + TRGT_STR2)/2) %>%
  mutate(Whole.length = RU.length*mean_repeat_count) %>% #head()
  group_by(STR,MOTIFS) %>%
  summarise(RU.length = mean(RU.length),mean_Whole.length = mean(Whole.length),mean_repeat_count= mean(mean_repeat_count)) %>% #head()
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) %>%
  ggplot(aes(x=RU.length,y=mean_repeat_count,color=GC)) +
  geom_point() 
  
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

  

#str_class_ref <- read_table("~/Desktop/KU/@research/STR/db/hg38_2020_rmsk.bed",col_names = F)
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
