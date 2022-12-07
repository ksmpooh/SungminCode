library(tidyverse)
library(stringr)
library(readxl)
library(ggbreak)
library(ggpubr)
library(dplyr)

setwd("~/Desktop/KCDC/HLA_seq/")

df <- read_xlsx("variant.xlsx",sheet = 2)
head(df)

colnames(df)
df %>% filter(TOOL == "GATK") %>% #head() 
  #filter(Seqeuncing )
  select(Seqeuncing,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`,...12) %>%  
  mutate(QC = ifelse(is.na(...12),QC,paste0(QC,"_",...12))) %>% #head()
  group_by(Seqeuncing) %>%
  pivot_longer(cols = c(RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%#  head()
  ggplot(aes(x=type,y=Val,color = QC)) +
  geom_point() +
  geom_line() + 
  #geom_bar() +
  facet_grid(~Seqeuncing)
  

df %>% filter((is.na(...12))) %>%
  select(Seqeuncing,TOOL,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  pivot_longer(cols = c(RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%  #head()
  ggplot(aes(x=type,y=Val,color = QC)) +
  geom_point() +
  geom_line() + 
  facet_grid(TOOL~Seqeuncing)


#aes(x = key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black'))))

df %>% filter((is.na(...12))) %>% #count(TOOL,QC)
  select(Seqeuncing,TOOL,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  pivot_longer(cols = c(RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%  #head()
  mutate(Method = paste0(TOOL,"_",QC)) %>% #head()
  mutate(Method = str_replace(Method, "Deepvariant-GLnexus", "DV")) %>% #count(Method)
  #ggplot(aes(x=Method,y=Val,fill = type)) +
  ggplot(aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),y=Val,fill = type)) +
  geom_bar(stat = "identity") + 
  theme(legend.title = element_blank()) +
  xlab(element_blank()) + ylab(element_blank()) + 
  scale_fill_discrete(guide="none") +
  facet_grid(~Seqeuncing) ->a1

  #geom_point() #+
  #geom_line() + 


df %>% filter((is.na(...12))) %>% #count(TOOL,QC)
  select(Seqeuncing,TOOL,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  pivot_longer(cols = c(RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%  #head()
  mutate(Method = paste0(TOOL,"_",QC)) %>% #head()
  mutate(Method = str_replace(Method, "Deepvariant-GLnexus", "DV")) %>% #head()
  ggplot(aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),y=Val,fill = type)) +
  geom_bar(position = 'fill',stat = "identity") + 
  theme(legend.title = element_blank(),legend.position="bottom") +
  xlab(element_blank()) + ylab(element_blank()) + 
  facet_grid(~Seqeuncing) ->a2


ggarrange(a1,a2,ncol = 1,nrow = 2,heights = c(8,9))
