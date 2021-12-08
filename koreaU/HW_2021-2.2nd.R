### 고려대학원 2021-2 : ㄷ데이터시각화 과제
library(tidyverse)
library(readxl)
library(grid)
library(gridExtra)
setwd("~/Desktop/KU/2021_Fall/biovisual/hw/ScienceDirect_files_19Oct2021_06-39-57.618/")

df1<- read_excel("1-s2.0-S0092867421008576-mmc1.xlsx",sheet = 2)

df1 %>% select(Participant,Type,ends_with("mutation")) %>% 
  #colnames()
  count(Participant)

df1 %>% select(ends_with("mutation")) %>% 
  count(TP53.mutation)

ref <- read_excel("1-s2.0-S0092867421008576-mmc1.xlsx",sheet = 3,skip = 1)
head(ref)
colnames(df)
df %>% select(Participant,Type,ends_with("mutation")) %>% 
  colnames()

df %>% select(ends_with("mutation")) %>% head()

  
  
  
df <- read_excel("1-s2.0-S0092867421008576-mmc2.xlsx",sheet = 10)
df %>% colnames()
a <-df %>% select(Sample.ID) %>% count(Sample.ID)
genelist <- df %>% select(Sample.ID,gene,type) %>%
  filter(grepl('C3L|C3N', Sample.ID)) %>% #count(Sample.ID)#head()
  filter(gene %in% ref$`Gene Symbol`)  %>% ##head()
  count(gene) %>% 
  filter(n > 2) %>% select(gene)


df %>% select(Sample.ID,gene,type) %>%
  filter(grepl('C3L|C3N', Sample.ID)) %>% #count(Sample.ID)#head()
  filter(gene %in% ref$`Gene Symbol`)  %>% ##head()
  filter(gene %in% genelist) %>%
  #pivot_wider(names_from = gene,values_from = type) %>% #head()
  ggplot(aes(x=Sample.ID,y=gene,fill=type)) +
    geom_tile()
  #spread(key = Sample.ID,value = gene) %>% head()


genelist <-c("KRAS","NF1","NRAS","HRAS","PIK3CA","PIL3R1",
            "AKT1","PTEN","NFE2L2","KEAP1","CUL3","ERBB2","EGFR",
            "FGFR2","FGFR3","FGFR1","RB1","CDKN2A","CCND1","NSD3","KMT2D",
            "KAT6A","ARID1A","SOX2","TP63","FAT1","TP53")

## gene : 27개
a <- df %>%  count(Gene) %>%
  filter(Gene %in% ref$`Gene Symbol`) %>%
  filter(n > 7)
a


## 차
genelist <-c("KRAS","NF1","NRAS","HRAS","PIK3CA","PIL3R1",
             "AKT1","PTEN","NFE2L2","KEAP1","CUL3","ERBB2","EGFR",
             "FGFR2","FGFR3","FGFR1","RB1","CDKN2A","CCND1","NSD3","KMT2D",
             "KAT6A","ARID1A","SOX2","TP63","FAT1","TP53")
library(tidyverse)
library(readxl)
library(grid)
library(gridExtra)
setwd("~/Desktop/KU/2021_Fall/biovisual/hw/ScienceDirect_files_19Oct2021_06-39-57.618/")

df7<- read_excel("1-s2.0-S0092867421008576-mmc7.xlsx",sheet = 7)
df7 %>% select(Variant_Type,Variant_Function) %>%
  #count(Variant_Type)
  count(Variant_Function)
  head()

df1<- read_excel("1-s2.0-S0092867421008576-mmc1.xlsx",sheet = 2)

df1 %>% select(Participant,Type,ends_with("mutation")) %>% 
  #colnames()
  count(Participant)

ref <- read_excel("1-s2.0-S0092867421008576-mmc1.xlsx",sheet = 3,skip = 1)
head(ref)
colnames(df)


df <- read_excel("1-s2.0-S0092867421008576-mmc2.xlsx",sheet = 10)
df %>% select(Sample.ID,gene,type,Mutation_Status,is_coding,is_flank,MUTATION_HOTSPOT,Variant_Classification) %>%
  #head()
  count(Variant_Classification)
  #count(MUTATION_HOTSPOT)
  #count(Mutation_Status)
  count(type)
  count(is_flank)
  head()

df %>% select(Sample.ID,gene,type) %>%
  filter(grepl('C3L|C3N', Sample.ID)) %>% #count(Sample.ID)#head()
  filter(gene %in% ref$`Gene Symbol`)  %>% ##head()
  filter(gene %in% genelist) %>%
  #pivot_wider(names_from = gene,values_from = type) %>% #head()
  ggplot(aes(x=Sample.ID,y=factor(gene,levels = rev(genelist)),fill=type)) +
  #ggplot(aes(x=Sample.ID,y=gene,fill=factor(type,levels = genelist))) +
  geom_tile(stat="identity") + 
  labs(y = "",x="") + 
  theme(legend.position ="bottom") + 
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

reve
  
 #spread(key = Sample.ID,value = gene) %>% head()



### 
ref2 <- read.table("")