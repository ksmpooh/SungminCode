library(stringr)
library(tidyverse)
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/michigan/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_han/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_pan/03.allele.matching/")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla/5M_28_33/03.allele.matching")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/03.allele.matching/")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_pan/8digit/03.allele.matching/")
#flist <- list.files("./",pattern = "missINFO.txt", invert=TRUE, value=TRUE)
#list.files("./",pattern = "missINFO.txt", invert=TRUE, value=TRUE)



################
flist = grep(list.files("./"),pattern = "missINFO.txt", invert=TRUE, value=TRUE)

flist

df <-read.table(flist[1],header = T)
a <- df %>% summarise(across(colnames(df)[-1],sum))
for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
  a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
}
a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum),"filename" = flist[1])
#str_split_fixed(str_split_fixed(flist[1],'\\.',5)[,3],"_",2)[,2]
out <- a

for (i in 2:length(flist)) {
  df <-read.table(flist[i],header = T)
  a <- df %>% summarise(across(colnames(df)[-1],sum))
  for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
    a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
  }
  a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum),"filename" = flist[i])
  out <- rbind(out,a)
}
head(out)

out %>% mutate(filename = str_replace_all(filename,".txt","")) %>%# head()
  mutate('digit' = ifelse(grepl(pattern = "td",filename),"2","4")) %>%
  mutate("CV" = str_split_fixed(filename,"\\.",3)[,2]) %>% #head()
  mutate("Ref" = str_split_fixed(filename,"\\.",4)[,4]) -> out

head(out)
out$Tool <- "Minimac4"
out$Tool <- "SNP2HLA"

#michigan <- out
#snp2hla1 <- out
snp2hla2 <- out
snp2hla_han <- out
snp2hla_pan <- out

snp2hla_pan1 <- out
snp2hla_pan2 <- out

head(michigan)

snp2hla %>% summarise(across(colnames(df)[-1],mean))

snp2hla %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref) %>% #count(Ref)
  filter(Ref %in% c("cmp_Nomencleaner.fdvstd","cmp_Nomencleaner")) %>% #head()
  filter(!(digit == "2" & Ref == "cmp_Nomencleaner")) %>% #head()
  mutate(Test = str_split_fixed(CV,"_",2)[,1]) %>%
  pivot_longer(1:9,names_to = 'Gene',values_to = 'Accuracy') %>%
  mutate(Gene = ifelse(grepl("HLA",Gene),str_split_fixed(Gene,"_",2)[,2],'Overall')) %>% #head()
  ggplot(aes(x=Gene,y=Accuracy,fill=Test)) +
  geom_bar(position='dodge', stat='identity') + 
  coord_cartesian(ylim = c(0.7, 1)) +
  facet_grid(~digit)
  
snp2hla %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref) %>% #count(Ref)
  filter(Ref %in% c("cmp_Nomencleaner.fdvstd","cmp_Nomencleaner")) %>% #head()
  filter(!(digit == "2" & Ref == "cmp_Nomencleaner")) %>% #head()
  mutate(Test = str_split_fixed(CV,"_",2)[,1]) %>%#-> a
  summarise(summarise(across(1:9,mean)))
  



head(michigan)
head(snp2hla)
out <- michigan
out <- rbind(michigan,snp2hla)
#out <- rbind(michigan,snp2hla,snp2hla_han,snp2hla_pan)

###############################

out %>% 
  filter(Ref %in% c("cmp_Nomencleaner.fdvstd","cmp_Nomencleaner")) %>%
  filter(!(digit == "2" & Ref == "cmp_Nomencleaner")) %>% #head()
  mutate(Test = str_split_fixed(CV,"_",2)[,1]) ->out
  

colnames(out)
out %>% #filter(digit == 2) %>%# head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref) %>%
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  #filter(digit == 2) %>%
  filter(digit == 4) %>%
  #group_by(CV) %>%
  ggplot(aes(x=Gene,y=Accuracy,color=CV,group=CV)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  facet_grid(~Ref)

out %>% #filter(digit == 2) %>%# head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref) %>%
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  filter(digit == 4) %>%
  ggplot(aes(x=Gene,y=Accuracy,color=CV)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  facet_grid(~Ref)

  
out %>% #filter(digit == 2) %>%# head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref) %>%
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  filter(digit == 2) %>%
  ggplot(aes(x=Gene,y=Accuracy,color=CV)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  facet_grid(~Ref)



out %>% #filter(digit == 2) %>%# head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref,Tool) %>%
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  filter(digit == 2) %>%
  ggplot(aes(y=Accuracy,x=Gene,fill=Tool,color=CV)) + 
  geom_bar(position='dodge', stat='identity') +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  facet_grid(~Ref)

out %>% #filter(digit == 2) %>%# head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref,Tool) %>%
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DRB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  filter(digit == 2) %>%
  ggplot(aes(y=Accuracy,x=Gene,color=CV)) + 
  geom_point() + 
  geom_line() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  facet_grid(~Ref)


###
head(michigan)
head(out)
table(out$Ref)
out %>% #filter(digit == 2) %>%# head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref,Tool) %>% #count(digit,Tool)#head()#count(CV)
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>% #count(Tool,Ref)#head()
  filter(Ref != "cmp_RealNGStyping") %>%  #count(Ref)
  filter(!(digit == 2 & Ref == "cmp_Nomencleaner")) %>% #count(Ref)
  mutate(Ref = ifelse(Tool == "Minimac4","Multi-ethnic","KMHC"),CV=str_split_fixed(CV,"_",2)[,1]) %>%
  mutate(Gene = ifelse(grepl("HLA",Gene),str_split_fixed(Gene,"_",2)[,2],'Overall')) %>% #head()
  filter(CV != '520sample') -> out1 # %>% #count(digit,CV,Ref)#head()#count(CV)

snp2hla_han$Ref <- "Han Chinese"
snp2hla_pan$Ref <- "Pan-Kor"

snp2hla_han %>% rbind(snp2hla_pan) %>% #head()#count(digit)
  filter(!grepl(filename,pattern = "_td")) %>% #count(digit)
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref,Tool) %>% #count(digit,Tool)#head()#count(CV)
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>% #count(Tool,Ref)#head()
  mutate(CV=str_split_fixed(CV,"_",2)[,1]) %>%
  mutate(Gene = ifelse(grepl("HLA",Gene),str_split_fixed(Gene,"_",2)[,2],'Overall')) %>% #head()
  filter(CV != '520sample') -> out2 # %>% #count(digit,CV,Ref)#head()#count(CV)

head(out1)
out1 %>% count(digit,Ref,Tool,CV)
head(out2)
out2 %>% count(digit,Ref,Tool)
out1 %>% count(Gene)
out2 %>% count(Gene)

  #count(filename)3
out1 %>% rbind(out2) %>% #head()#count(CV)
  filter(!(Gene %in% c("A","B","DRB1"))) %>%
  filter(Gene != "Overall") %>%
  #filter(Gene %in% c("A","B","DRB1")) %>%
  mutate(digit = paste0(digit,' digit')) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=factor(Ref,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")))) +
  geom_boxplot() + 
  theme(legend.title = element_blank()) +
  labs(fill = "Reference panel") + 
  theme(legend.position = "bottom") + 
  theme(legend.title=element_text(size=13,face = "bold"),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14,face = "bold"))+
  facet_grid(~digit) +
  theme(strip.text.x = element_text(size = 13,face = "bold"))

out1 %>% rbind(out2) %>% #head()#count(CV)
  filter(Gene == "Overall") %>% #head()
  group_by(Ref,digit) %>%
  summarise(Accuracy = mean(Accuracy)*100) %>%
  mutate(digit = paste0(digit,' digit')) %>% #head()
  #arrange(Ref  = factor(Ref,levels = c("KMHC","Multi-Ethnic","Han Chinese","Pan-Kor"))) %>%
  ggplot(aes(x=factor(Ref,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")),y=Accuracy,fill=factor(Ref,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")))) +
  geom_bar(position="dodge",stat = "identity") + 
  geom_text(size = 4,aes(label = paste(round(Accuracy,2),"%")),hjust = 0.5, vjust = 2, position = "stack") +
  theme(legend.title = element_blank()) +
  labs(fill = "Reference panel") + 
  ylab("Accuracy (%)") +
  theme(legend.position = "bottom") + 
  theme(legend.title=element_text(size=13,face = "bold"),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14,face = "bold"))+
  facet_grid(~digit) +
  theme(strip.text.x = element_text(size = 13,face = "bold"))


out1 %>% rbind(out2) %>% #head()
  pivot_wider(names_from = Gene,values_from = Accuracy) %>%
  select(-CV,-Tool) %>%
  group_by(digit,Ref) %>% #head()
  summarise(across(1:9,mean)) -> a

head(snp2hla_pan)
head(snp2hla_han)

snp2hla_han %>% rbind(snp2hla_pan) %>% #head()#count(digit)
  filter(!grepl(filename,pattern = "_td")) %>% #count(digit)
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,digit,CV,Ref,Tool) %>% #count(digit,Tool)#head()#count(CV)
  pivot_longer(cols = c(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall),names_to = 'Gene',values_to = 'Accuracy') %>% #count(Tool,Ref)#head()
  mutate(CV=str_split_fixed(CV,"_",2)[,1]) %>%
  mutate(Gene = ifelse(grepl("HLA",Gene),str_split_fixed(Gene,"_",2)[,2],'Overall')) %>% #head()
  #filter(CV == '520sample') %>% 
  #pivot_wider(names_from = Gene,values_from = Accuracy) -> a
  filter(Ref == "Pan-Kor") %>%
  ggplot(aes(x= Gene,y=Accuracy,fill=CV)) +
  geom_bar(position='dodge', stat='identity') +
  facet_grid(~digit)
  

###### 2 digit : only 2 vs 2 for 4
head(michigan)
head(snp2hla)
#head(snp2hla_han)
#head(snp2hla_pan)

michigan %>% rbind(snp2hla) %>% filter(digit == 2) %>%
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,Tool) %>% #count(CV,Tool,Ref)#head()#count(CV)
  filter(Ref != "cmp_RealNGStyping") %>% 
  mutate('from' = ifelse(Ref == "cmp_Nomencleaner","2digit","4to2digit")) %>%
  mutate(Tool=ifelse(Tool=="Minimac4","Michigan(Multi-ethnic) : Minimac4","HLA-TAPAS(KMHC) : SNP2HLA")) %>%
  filter(CV != "520sample") %>%
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=from))+
  geom_boxplot()+
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~Tool, ncol = 1) +  
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))
  

###### KMHC 10M vs 5M
head(snp2hla1) #5M
head(snp2hla2) #10M

snp2hla1$region <- "10M"
snp2hla2$region <- "5M"

snp2hla1 %>% rbind(snp2hla2) %>%
  filter(grepl("cmp_Nomencleaner",Ref)) %>% 
  filter(!(digit == "2" & Ref == "cmp_Nomencleaner")) %>% #head()
  mutate(CV = str_replace_all(CV,"_5CV","")) %>%#count(CV) #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Tool,region,digit) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = "Accuracy") %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  mutate(digit = paste0(digit,' digit')) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=region))+
  geom_boxplot()+
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~digit, ncol = 1) +  
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))


#### KMHC IMGT 3220 vs 3370

imgt3320 <- out
imgt3320$IMGT <- "3320"

imgt3370 <- out
imgt3370$IMGT <- "3370"

head(imgt3320)
head(imgt3370)

out <- rbind(imgt3320,imgt3370)

head(out)


out %>% filter(digit == 2) %>% #head()
  select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,IMGT) %>% #count(CV,Tool,Ref)#head()#count(CV)
  filter(Ref != "cmp_RealNGStyping") %>% 
  mutate('from' = ifelse(Ref == "cmp_Nomencleaner","2digit","4to2digit")) %>%
  #mutate(Tool=ifelse(Tool=="Minimac4","Michigan(Multi-ethnic) : Minimac4","HLA-TAPAS(KMHC) : SNP2HLA")) %>%
  filter(CV != "520sample") %>%
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=from))+
  geom_boxplot() +
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~IMGT, ncol = 2) +  
  theme(legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))

#out %>% filter(IMGT == "3320") %>% #head()#count()#count(digit)
out %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,IMGT) %>% #count(CV,Tool,Ref)#head()#count(CV)
  filter(Ref != "cmp_RealNGStyping") %>% #head()#count(digit) 
  filter(!(digit == 2 & Ref != "cmp_Nomencleaner.fdvstd")) %>% #head()#count()
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=IMGT))+
  geom_boxplot() +
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~digit, ncol = 2) +  
  theme(#legend.title=element_blank(),
        legend.text=element_text(size=11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))


out %>% head()
out %>% select(HLA_A.empty,HLA_B.empty,HLA_C.empty,HLA_DRB1.empty,HLA_DPA1.empty,HLA_DPB1.empty,HLA_DQA1.empty,HLA_DQB1.empty,empty_Sum,CV,Ref,digit,IMGT) %>% #count(CV,Tool,Ref)#head()#count(CV)
  filter(Ref != "cmp_RealNGStyping") %>% #head()#count(digit) 
  filter(!(digit == 2 & Ref != "cmp_Nomencleaner.fdvstd")) %>% #head()#count()
  pivot_longer(1:9,names_to = "Gene",values_to = 'Count') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  mutate(Gene = str_replace_all(Gene,".empty","")) %>% #head()
  mutate(Gene = ifelse(Gene == "empty_Sum","overall",Gene)) %>%
  group_by(IMGT,digit,Gene) %>% #head()
  #count(Count)
  summarise("NA_count_sum" = sum(Count)) %>%
  ggplot(aes(x=Gene,y=NA_count_sum,color=IMGT))+
  #geom_boxplot() +
  geom_point() +
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~digit, ncol = 2) +  
  theme(#legend.title=element_blank(),
    legend.text=element_text(size=11),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))



head(out)

df <- read
  
  