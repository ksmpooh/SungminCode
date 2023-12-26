library(stringr)
library(tidyverse)

setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/michigan/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_han/03.allele.matching/")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_pan/03.allele.matching/")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla/5M_28_33/03.allele.matching")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/03.allele.matching/")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320_2to4/03.allele.matching/")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320_2to4_5M/03.allele.matching/")
setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/TEST.JG/03.allele.matching/")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_1kgp/03.allele.matching/")
#setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_pan/8digit/03.allele.matching/")
#flist <- list.files("./",pattern = "missINFO.txt", invert=TRUE, value=TRUE)
#list.files("./",pattern = "missINFO.txt", invert=TRUE, value=TRUE)

setwd("/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320/03.allele.matching_freq/")

################
flist = grep(list.files("./"),pattern = "missINFO", invert=TRUE, value=TRUE)

flist
head(df)
head(a)
df <-read.table(flist[1],header = T)
head(df)
a <- df %>% summarise(across(colnames(df)[-1],sum))
#for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
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
#out$Tool <- "SNP2HLA"

michigan <- out
#snp2hla1 <- out
#snp2hla <- out
#snp2hla2 <- out
#snp2hla3 <- out
#snp2hla2 <- out
#snp2hla_han <- out
snp2hla_pan <- out

snp2hla_pan1 <- out
snp2hla_pan2 <- out

head(michigan)




snp2hla <- out
snp2hla <- out
snp2hla <- out

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
  

head(snp2hla_pan)

head(michigan)
head(snp2hla)
out <- michigan
out <- rbind(michigan,snp2hla)
out <- rbind(snp2hla,snp2hla2)
out <- rbind(snp2hla,snp2hla2,snp2hla3)
out <- rbind(michigan,snp2hla,snp2hla_han,snp2hla_pan)
head(snp2hla)
head(snp2hla2)
###############################
head(out)
out %>% 
  filter(Ref %in% c("cmp_Nomencleaner.fdvstd","cmp_Nomencleaner")) %>% #head()
  filter(!(digit == "2" & Ref == "cmp_Nomencleaner")) %>% #count(Tool,Ref)
  mutate(Test = str_split_fixed(CV,"_",2)[,1]) ->out
  
head(out)
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
  #filter(!(Gene %in% c("A","B","DRB1"))) %>%
  #filter(Gene != "Overall") %>%
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
  #facet_grid(~digit) +
  facet_wrap(~digit,ncol = 1)
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
  

head(out)
out %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,filename,digit,CV,Ref,Tool,Test) %>%
  mutate(panel = str_split_fixed(filename,"\\.",5)[,3]) %>% 
  mutate(panel = ifelse(panel == "michiganHLAimp_fd","Multi-ethnic",ifelse(panel=="SNP2HLAHLAimp_fd","KMHC",ifelse(panel=="SNP2HLAHLAimp_HanREF_fd","Han Chinese","Pan-Kor")))) %>%
  filter(Test != "520sample") %>% 
  select(-filename,-Ref,-CV,-Tool) %>% 
  pivot_longer(1:8) %>% #head()
  group_by(digit,panel,name) %>% 
  mutate(name = str_split_fixed(name,"_",2)[,2]) %>%
  filter(!name %in% c("A","B","DRB1")) %>%
  summarise(mean = mean(value)) %>%
  ggplot(aes(x=name,y=mean,fill=factor(panel,levels = c("KMHC","Multi-ethnic","Han Chinese","Pan-Kor")))) +
  geom_bar(position="dodge",stat = "identity") + 
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

##### kMHC ori vs HLA type 2for4

head(snp2hla)
head(snp2hla2)

snp2hla$type = "IMGT3320_original"
snp2hla2$type = "IMGT3320_modify"


snp2hla %>% rbind(snp2hla2) %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,type) %>% #count(CV,Tool,Ref)#head()#count(CV)
  filter(Ref != "cmp_RealNGStyping") %>% #head()#count(digit) 
  filter(!(digit == 2 & Ref != "cmp_Nomencleaner.fdvstd")) %>% #head()#count()
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=type))+
  geom_boxplot() +
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~digit, ncol = 2) +  
  theme(#legend.title=element_blank(),
    legend.text=element_text(size=11),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))


df <-read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/snp2hla_imgt3320_2to4/01.impresult/test.txt")
head(f)
  

#####
head(snp2hla)
head(snp2hla2)

snp2hla$type = "IMGT3320_original"
snp2hla2$type = "IMGT3320_modify_10M"
snp2hla3$type = "IMGT3320_modify_5M"


snp2hla %>% rbind(snp2hla2) %>% rbind(snp2hla3) %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,Ref,digit,type) %>% #count(CV,Tool,Ref)#head()#count(CV)
  filter(Ref != "cmp_RealNGStyping") %>% #head()#count(digit) 
  filter(!(digit == 2 & Ref != "cmp_Nomencleaner.fdvstd")) %>% #head()#count()
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Accuracy,fill=type))+
  geom_boxplot() +
  #facet_grid(~Tool,rows = vars(Tool))
  facet_wrap(~digit, ncol = 2) +  
  theme(#legend.title=element_blank(),
    legend.text=element_text(size=11),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14,face = "bold"))+
  theme(strip.text.x = element_text(size = 13,face = "bold"))


#### VS JG 2digit

flist = grep(list.files("./"),pattern = "missINFO.txt", invert=TRUE, value=TRUE)
flist
df <-read.table(flist[1],header = T)
a <- df %>% summarise(across(colnames(df)[-1],sum))
for (gene in c("HLA_A","HLA_B","HLA_DRB1")) {
  a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
}
a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum),"filename" = flist[1])
#str_split_fixed(str_split_fixed(flist[1],'\\.',5)[,3],"_",2)[,2]
out <- a
head(out)
for (i in 2:length(flist)) {
  df <-read.table(flist[i],header = T)
  a <- df %>% summarise(across(colnames(df)[-1],sum))
  for (gene in c("HLA_A","HLA_B","HLA_DRB1")) {
    a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
  }
  a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum),"filename" = flist[i])
  out <- merge(out,a,all=T)
}
head(out)
a <- out

out <- a
out %>% mutate(filename = str_replace_all(filename,".txt","")) %>%# head()
  mutate('digit' = ifelse(grepl(pattern = "fd",filename),"4to2","2")) %>%
  mutate("CV" = str_split_fixed(filename,"\\.",3)[,2]) %>% #head()
  #mutate("digit" = ifelse(grepl("fd"str_split_fixed(filename,"\\.",8)[,8]))),"2","4to2")  %>%
  mutate("Ref" = str_split_fixed(filename,"\\.",8)[,5]) -> out

head(out)


out %>% select(HLA_A,HLA_B,HLA_DRB1,overall,Ref,digit) %>% #head()
  filter(digit == "4to2") %>%
  pivot_longer(1:4) %>% mutate(Ref=str_split_fixed(Ref,"_",2)[,1]) %>% #head()
  ggplot(aes(x=name,y=value,color=Ref)) + 
  geom_point()
  #facet_grid(~digit)


###### 20231121 by frequency

head(out)
out$freq <- str_split_fixed(out$Ref,"\\.",2)[,2]
out %>% count(freq)

out %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,freq) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% #head()
  filter(freq != "notrare") %>% #count(freq,Gene)
  group_by(freq,Gene) %>% #head()
  summarise(Accuracy = round(mean(Accuracy)*100,1)) -> a
  
head(out)
#writexl::write_xlsx(out,"~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.5CV.result.byfreq.xlsx")



out %>% select(HLA_A,HLA_B,HLA_C,HLA_DRB1,HLA_DPA1,HLA_DPB1,HLA_DQA1,HLA_DQB1,overall,CV,freq) %>% #count(CV,Tool,Ref)#head()#count(CV)
  pivot_longer(1:9,names_to = "Gene",values_to = 'Accuracy') %>% #head()
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% #head()
  #filter(freq != "notrare") %>% #filter(is.na(Accuracy))
  na.omit() %>%
  group_by(freq,Gene) %>% #head()
  summarise(Accuracy = round(mean(Accuracy)*100,1)) %>% #head()
  ungroup() %>%
  ggplot(aes(x=fct_rev(freq),y=Accuracy,fill=freq))+
  geom_bar(stat='identity',position = 'dodge')  + 
  facet_wrap(~Gene, ncol = 3) +  
  geom_text(size = 4,aes(label = Accuracy),hjust = 1, vjust = 0.5,position = position_dodge(width = .9)) +
  scale_fill_discrete(name="Frequency") +
  theme(#legend.title=element_text("bold"),
      #legend.title=element_text("Ref.panel"),
      legend.text=element_text(size=11),
      legend.position = "bottom",
      axis.text.y = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()) + 
    theme(strip.text.x = element_text(size = 13,face = "bold")) +
    coord_flip(clip = "on",ylim = c(0, 100)) #-> p4
  


df <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Result/KMHC.5CV.result.byfreq.xlsx")

colnames(df)
df %>% select(grep(".match",colnames(df)),grep(".wrong",colnames(df)),match_Sum,wrong_Sum,CV,freq) %>% #head()
  rename('all.match' = match_Sum, 'all.wrong' = wrong_Sum) %>% #head()
  group_by(freq) %>% #colnames()
  #summarise(across(HLA_A.match:HLA_DQB1.wrong, sum)) %>% #head()
  summarise(across(HLA_A.match:all.wrong, sum)) -> df

head(df)

df %>% select(freq,grep(".match",colnames(df))) %>%
  pivot_longer(2:10,names_to = "Gene",values_to = 'match_count') %>% 
  mutate(Gene = str_split_fixed(Gene,"\\.",2)[,1]) -> a

df %>% select(freq,grep(".wrong",colnames(df))) %>%
  pivot_longer(2:10,names_to = "Gene",values_to = 'wrong_count') %>% 
  mutate(Gene = str_split_fixed(Gene,"\\.",2)[,1]) -> b

a %>% left_join(b) %>% filter(freq %in% c("common","lesscommon")) %>% #head()
  group_by(Gene) %>%
  summarise(across(match_count:wrong_count,sum)) %>% mutate(freq = "not_Rare")  -> c
head(a)  
head(c)

a %>% left_join(b) %>% rbind(c) %>% mutate(Accuracy = match_count/(match_count+wrong_count), Allele_count = paste0(match_count,"/",match_count + wrong_count)) %>% #head()
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% #head()
  mutate(Gene = ifelse(Gene == "all","Overall",Gene)) %>% #head()
  filter(freq != "notrare") %>% #head()
  na.omit() %>% #head()
  mutate(Accuracy = round((Accuracy)*100,1)) %>%
  ggplot(aes(x=fct_rev(freq),y=Accuracy,fill=freq))+
  geom_bar(stat='identity',position = 'dodge')  + 
  facet_wrap(~Gene, ncol = 3) +  
  geom_text(size = 5,aes(label = Accuracy),hjust = 1, vjust = 0.5,position = position_dodge(width = .9)) +
  geom_text(size = 4,aes(label = Allele_count),y=00,hjust = 0) +
  scale_fill_discrete(name="Frequency") +
  theme(#legend.title=element_text("bold"),
    #legend.title=element_text("Ref.panel"),
    legend.text=element_text(size=11),
    legend.position = "bottom",
    axis.text.y = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()) + 
  theme(strip.text.x = element_text(size = 13,face = "bold")) +
  #coord_flip(clip = "on",ylim = c(0, 100)) #-> p4
  coord_flip(clip = "on")

