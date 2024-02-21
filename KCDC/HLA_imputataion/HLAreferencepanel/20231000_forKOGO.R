library(stringr)
library(tidyverse)
library(ggpubr)
#setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/rmnotincludetype/maf/")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/rmnotincludetype/v1")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/rmnotincludetype/v2")

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320_maf/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/rmnotincludetype/maf/")


flist_all = grep(list.files("./"),pattern = "missINFO.txt", invert=TRUE, value=TRUE)
flist_all
#flist = flist_all[grep(flist_all,pattern="td")]
flist = flist_all[grep(flist_all,pattern="fd")]
flist
#df %>% summarise(across(colnames(df)[-1],sum))
out <- NULL
for (i in 1:length(flist)) {
  df <-read.table(flist[i],header = T)
  a <- df %>% summarise(across(colnames(df)[-1],sum))
  for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
    a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
  }
  a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum),"filename" = flist[i])
  out <- rbind(out,a)
}
head(out)
dim(out)
out %>% mutate(filename = str_replace_all(filename,".txt","")) %>%# head()
  mutate('digit' = ifelse(grepl(pattern = "td",filename),"2","4")) %>%
  mutate("Ref" = str_split_fixed(filename,"\\.",4)[,4]) %>%
  mutate("panel" = str_split_fixed(filename,"\\.",4)[,3])  %>%
  mutate("CV" = str_split_fixed(filename,"\\.",3)[,2]) %>% #head()
  mutate(panel = str_replace_all(panel,"_fd","")) %>% #head()
  mutate(panel = str_replace_all(panel,"_td","")) -> out
head(out)
table(out$Ref)
table(out$CV)
'''
out %>% 
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp", "SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF", "michiganHLAimp"), labels = c("KMHC", "Han Chinese", "Pan-Kor", "Multi-ethnic"))) %>% 
  writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/HLAimputation_result_compare.4panels.xlsx")
'''


out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  #filter(Ref != "cmp_RealNGStyping") %>% 
  #mutate('from' = ifelse(Ref == "cmp_Nomencleaner","2digit","4to2digit")) %>%
  #mutate(Tool=ifelse(Tool=="Minimac4","Michigan(Multi-ethnic) : Minimac4","HLA-TAPAS(KMHC) : SNP2HLA")) %>%
  #filter(CV != "520sample") %>%
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #count(panel)
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp", "SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF", "michiganHLAimp"), labels = c("KMHC", "Han Chinese", "Pan-Kor", "Multi-ethnic"))) %>% 
  ggplot(aes(x=Gene,y=Accuracy,fill=panel))+
  geom_boxplot() +
  facet_wrap(~digit, ncol = 2) +  
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))
head(out)

out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  group_by(panel,digit) %>% #head()
  summarise(Accuracy = mean(Accuracy)) %>%
  #mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(digit = ifelse(digit == "2","One-field","Two-field")) %>% 
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp","michiganHLAimp","SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  ceiling(Accuracy)
  

out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  group_by(panel,digit) %>% #head()
  summarise(Accuracy = mean(Accuracy)) %>%
  #mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(digit = ifelse(digit == "2","One-field","Two-field")) %>%
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp","michiganHLAimp","SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  ggplot(aes(x=panel,y=Accuracy,fill=panel))+
  #geom_boxplot() +
  geom_bar(stat='identity') + 
  coord_cartesian(ylim = c(0.5, 1)) + 
  facet_wrap(~digit, ncol = 2) +  
  #geom_text(size = 4,aes(label = paste(round(Accuracy,2),"%")),hjust = 0.5, vjust = 2, position = "stack") +
  geom_text(size = 4,aes(label = paste(round(Accuracy,2),"%")),hjust = 0.5, vjust = 2, position = "stack") +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold")) #-> p1
head(out)

out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  group_by(panel,digit) %>% #head()
  #summarise(Accuracy = ceiling(mean(Accuracy)*100)) %>% #head()
  summarise(Accuracy = round(mean(Accuracy)*100,1)) %>% #head()
  #mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(digit = ifelse(digit == "2","One-field","Two-field")) %>%
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp","michiganHLAimp","SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  ungroup() %>% 
  ggplot(aes(x=panel,y=Accuracy,fill=panel))+
  #geom_boxplot() +
  geom_bar(stat='identity',position = 'dodge') + 
  facet_wrap(~digit, ncol = 2) +  
  geom_text(size = 4,aes(label = Accuracy),hjust = 1, vjust = 0.5,position = position_dodge(width = .9)) +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) + 
  theme(strip.text.x = element_text(size = 13,face = "bold")) +
  coord_flip(clip = "on",ylim = c(50, 100))-> p1

p1

out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  #summarise(Accuracy = mean(Accuracy)) %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(digit = ifelse(digit == "2","One-field","Two-field")) %>%
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp","michiganHLAimp","SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  mutate(class = ifelse(Gene %in% c("A","B","C"),"class I","class II")) %>% #head()
  group_by(panel,digit,class) %>%
  #summarise(Accuracy = round(mean(Accuracy)*100)) %>% #head()
  summarise(Accuracy = round(mean(Accuracy)*100,1)) %>% #head()
  #summarise(Accuracy = mean(Accuracy)) %>% #head()
  ungroup() %>%
  ggplot(aes(x=fct_rev(class),y=Accuracy,fill=panel))+
  #geom_boxplot() +
  geom_bar(stat='identity',position = 'dodge') + 
  facet_wrap(~digit, ncol = 2) +  
  geom_text(size = 4,aes(label = Accuracy),hjust = 1, vjust = 0.5,position = position_dodge(width = .9)) +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        legend.position = "none",
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) + 
  theme(strip.text.x = element_text(size = 13,face = "bold")) +
  coord_flip(clip = "on",ylim = c(40, 100))-> p2

out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(digit = ifelse(digit == "2","One-field","Two-field")) %>%
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp","michiganHLAimp","SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  filter(Gene %in% c("A","B","DRB1")) %>%
  group_by(panel,digit,Gene) %>% #head()
  #summarise(Accuracy = mean(Accuracy)) %>% #head()
  #summarise(Accuracy = round(mean(Accuracy)*100)) %>% #head()
  summarise(Accuracy = round(mean(Accuracy)*100,1)) %>% #head()
  ungroup() %>%
  ggplot(aes(x=fct_rev(Gene),y=Accuracy,fill=panel))+
  #geom_boxplot() +
  geom_bar(stat='identity',position = 'dodge') + 
  facet_wrap(~digit, ncol = 2) +  
  geom_text(size = 4,aes(label = Accuracy),hjust = 1, vjust = 0.5,position = position_dodge(width = .9)) +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        legend.position = "none",
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) + 
  theme(strip.text.x = element_text(size = 13,face = "bold")) +
  coord_flip(clip = "on",ylim = c(80, 100))-> p3
p3

out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% 
  mutate(digit = ifelse(digit == "2","One-field","Two-field")) %>%
  mutate(panel = factor(panel, levels = c("SNP2HLAHLAimp","michiganHLAimp","SNP2HLAHLAimp_HanREF", "SNP2HLAHLAimp_PanREF"), labels = c("KMHC", "Multi-ethnic","Han Chinese", "Pan-Kor"))) %>% 
  filter(Gene %in% c("A","B","DRB1")) %>%
  group_by(panel,digit,Gene) %>% #head()
  #summarise(Accuracy = mean(Accuracy)) %>% #head()
  #summarise(Accuracy = round(mean(Accuracy)*100)) %>% #head()
  summarise(Accuracy = round(mean(Accuracy)*100,1)) %>% #head()
  ungroup() %>%
  ggplot(aes(x=fct_rev(Gene),y=Accuracy,fill=panel))+
  #geom_boxplot() +
  geom_bar(stat='identity',position = 'dodge') + 
  facet_wrap(~digit, ncol = 2) +  
  geom_text(size = 4,aes(label = Accuracy),hjust = 1, vjust = 0.5,position = position_dodge(width = .9)) +
  scale_fill_discrete(name="Ref.Panel") +
  theme(#legend.title=element_blank(),
        #legend.title=element_text("Ref.panel"),
        legend.text=element_text(size=11),
        legend.position = "bottom",
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) + 
  theme(strip.text.x = element_text(size = 13,face = "bold")) +
  coord_flip(clip = "on",ylim = c(80, 100))-> p4

p4<- cowplot::get_legend(p4)
p4

p0 <- ggarrange(p1, p2, p3,widths = c(1,1.2,1.2),labels = c("A","B","C"),ncol = 3)
p0 <- annotate_figure(p0,
                #top = text_grob("Visualizing len", color = "red", face = "bold", size = 14),
                bottom = text_grob("Accuracy (%)", size = 15,face = "bold")
)

p <- ggarrange(p0, p4,heights = c(4,1),ncol = 1)
p
#annotate_figure(p,
                


gga            
             
head(out)

geom_text(aes(variable, `(all)`, label = sprintf("%2.1f", `(all)`), group = ustanova), 
          position = position_dodge(width = .9)) +



out %>% #filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_DRB1,CV,Ref,digit) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:3,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=Ref))+
  geom_boxplot() +
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~digit, ncol = 1) +  
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))


out %>% #

head(out)
out$QC <- "notpanel"
out$QC <- "onpanel"

out$QC <- "ori"
out$QC <- "ori_v1"
out$QC <- "ori_v2"
out$QC <- "maf"
out$QC <- "maf_v2"

#a <- out
#b <- out
#c <- out
#d <- out
e <- out
a %>% filter(panel == "SNP2HLAHLAimp") -> d
b %>% filter(panel == "SNP2HLAHLAimp") -> c
head(a)
head(b)


d %>% group_by(digit) %>%
  summarise(mean = mean(overall))

a %>% filter(panel == "SNP2HLAHLAimp") %>% 
  group_by(digit) %>%
  summarise(mean = mean(overall))


a %>% filter(panel == "SNP2HLAHLAimp") %>% 
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% #head()
  #mutate(Class = ifelse())
  #filter(Gene %in% c("A","B","C")) %>%
  #filter(Gene %in% c("DRB1","DPA1","DPB1","DQA1","DQB1")) %>%
  #filter(Gene %in% c("A","B","DRB1")) %>%
  group_by(digit,Gene) %>%
  summarise(mean = mean(Accuracy)) %>%
  pivot_wider(values_from = mean,names_from = Gene) #%>% writexl::write_xlsx("KMHC_accuracy.xlsx")

kmhc_freq1 %>% group_by() %>%
  mutate(common = ifelse(KMHC >= 0.05,"common",ifelse(KMHC < 0.01,"rare","less_common"))) %>%
  count(common) -> t


head(a)
head(b)

a %>% rbind(b) %>% rbind(c) %>% rbind(d) %>% rbind(e) %>% #head()
  filter(panel == "SNP2HLAHLAimp") %>%
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,panel,QC) %>% #count(CV,Tool,Ref)#head()#count(CV)
  #filter(Ref != "cmp_RealNGStyping") %>% 
  #mutate('from' = ifelse(Ref == "cmp_Nomencleaner","2digit","4to2digit")) %>%
  #mutate(Tool=ifelse(Tool=="Minimac4","Michigan(Multi-ethnic) : Minimac4","HLA-TAPAS(KMHC) : SNP2HLA")) %>%
  #filter(CV != "520sample") %>%
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=QC))+
  geom_boxplot() +
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~digit, nrow = 1) +  
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))


grep("empty",colnames(a))
a %>% rbind(b) %>% 
  filter(panel == "SNP2HLAHLAimp") %>% #head()
  select(colnames(a)[grep("wrong",colnames(a))],QC,CV,panel,digit) %>% #count(panel,QC,CV,digit)#head()
  pivot_longer(1:9) %>% #count(panel,QC,CV)#head()
  pivot_wider(names_from = QC,values_from = value) %>% #select(notpanel)
  ggplot(aes(x=notpanel,y=onpanel))+
  geom_point() + 
  facet_grid(~digit)

#hlatype <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLAtype_check/")  

head(a)


