#### QC ÈÄ ±×¸²


load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.PASS.intersect.RData")
load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.upper0.8.intersect.RData")

ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)
common_str <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/TRID_common_EHpass_TRGTupper0.8.txt",header = T)

head(trgt_pass_common)
head(eh_pass_common)


#### motif freqeuncy
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



#### EH vs TRGT length

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
            conf.int = TRUE, # Add confidence interval
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

