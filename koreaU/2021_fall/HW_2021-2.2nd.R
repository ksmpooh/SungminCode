### 고려대학원 2021-2 : ㄷ데이터시각화 과제
library(tidyverse)
library(readxl)
setwd("~/Desktop/KU/2021_Fall/biovisual/hw/ScienceDirect_files_19Oct2021_06-39-57.618/")

df1<- read_excel("1-s2.0-S0092867421008576-mmc1.xlsx",sheet = 2)
colnames(df1)
df1 %>% select()
df1 %>% select(Participant,Type,ends_with(".cna.wes"),ends_with('.pathway.alteration'),FGFR3.TACC3.fusion.rna) %>% 
#CDKN2A.pathway.alteration"
#PIK3CA.pathway.alteration" 
#"FGFR3.TACC3.fusion.rna"
  #colnames()
  mutate()
#  head()
#  count(Participant)
#df1 %>% select(Participant,ends_with('.pathway.alteration')) %>% 
a <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
#df1 %>% select(Participant,ends_with('.pathway.alteration')) %>% 
  #mutate(test1 = strsplit(CDKN2A.pathway.alteration,"\\|"))
  mutate(test1_1 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,1],
         test1_2 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,2],
         test2_1 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,1],
         test2_2 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,2],
         test2_3 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,3]) %>%
  #mutate(test = str_detect(test1_1,"amp$"))
  #mutate(test=)
  filter(str_detect(test1_1,"amp")|str_detect(test1_2,"amp")|str_detect(test2_1,"amp")|str_detect(test2_2,"amp")|str_detect(test2_3,"amp"))
  mutate(test = )

    str_extract_all(test1_1,pattern = regex(str_c("amp")),)
  
  #select(Participant,starts_with('test')) %>%
  #filter()
  extract(test1_1,c("Amplification","Deleltion","Mutation"),
          "%_amp")
  #head()
#  pivot_wider()

a <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>% 
   

  #str_subset(CDKN2A.pathway.alteration,"")
  
head()
  
  
  
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
test <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  #df1 %>% select(Participant,ends_with('.pathway.alteration')) %>% 
  #mutate(test1 = strsplit(CDKN2A.pathway.alteration,"\\|"))
  mutate(test1_1 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,1],
         test1_2 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,2],
         test2_1 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,1],
         test2_2 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,2],
         test2_3 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,3]) %>%
  #mutate(test = str_detect(test1_1,"amp$"))
  #mutate(test=)
  #filter(str_detect(test1_1,"amp")|str_detect(test1_2,"amp")|str_detect(test2_1,"amp")|str_detect(test2_2,"amp")|str_detect(test2_3,"amp"))
  filter(str_detect(test1_1,"amp")) %>%
  select(Participant,test1_1) %>%
  #str_replace_all(test1_1,"_amp$","")
  mutate(amplification = str_replace_all(test1_1,"_amp",""),type = "Amplication") %>% 
  pivot_wider(names_from = amplification,values_from = type)

test




test <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  #df1 %>% select(Participant,ends_with('.pathway.alteration')) %>% 
  #mutate(test1 = strsplit(CDKN2A.pathway.alteration,"\\|"))
  mutate(test1_1 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,1],
         test1_2 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,2],
         test2_1 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,1],
         test2_2 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,2],
         test2_3 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,3]) %>%
  #mutate(test = str_detect(test1_1,"amp$"))
  #mutate(test=)
 # filter(str_detect(test1_1,"amp")|str_detect(test1_2,"amp")|str_detect(test2_1,"amp")|str_detect(test2_2,"amp")|str_detect(test2_3,"amp")) %>%
#  select(Participant,starts_with("test"))
  filter(str_detect(test1_1,"amp")) %>%
  select(Participant,test1_1) %>%
  #str_replace_all(test1_1,"_amp$","")
  mutate(amplification = str_replace_all(test1_1,"_amp",""),type = "Amplication") %>% 
  pivot_wider(names_from = amplification,values_from = type)

test

make_new_df <- function(.data,df,inpattern,intype){
  .data %>% filter(str_detect(df,inpattern)) %>%
  select(Participant,df) %>%
  mutate(tmp = str_replace_all(df,paste0("_",inpattern),""),type = intype) %>% 
  select(Participant,tmp,type) %>%
  pivot_wider(names_from = tmp,values_from = type)
}  

test <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  mutate(test1_1 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,1],
         test1_2 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,2]) %>%
  make_new_df(test1_1,"amp","amplification")



test11 <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  mutate(test1_1 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,1],
         test1_2 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,2]) %>%
  filter(str_detect(test1_1,"amp")) %>%
  select(Participant,test1_1) %>%
  mutate(amplification = str_replace_all(test1_1,"_amp",""),type = "Amplication") %>% 
  select(Participant,amplification,type) %>%
  pivot_wider(names_from = amplification,values_from = type)

test12 <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  mutate(test1_1 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,1],
         test1_2 = str_split_fixed(CDKN2A.pathway.alteration,"\\|",3)[,2]) %>%
  filter(str_detect(test1_2,"amp")) %>%
  select(Participant,test1_2) %>%
  mutate(amplification = str_replace_all(test1_2,"_amp",""),type = "Amplication") %>% 
  select(Participant,amplification,type) %>%
  pivot_wider(names_from = amplification,values_from = type)


test21 <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  mutate(test2_1 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,1],
         test2_2 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,2],
         test2_3 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,3]) %>%
  filter(str_detect(test2_1,"amp")) %>%
  select(Participant,test2_1) %>%
  mutate(amplification = str_replace_all(test2_1,"_amp",""),type = "Amplication") %>% 
  select(Participant,amplification,type) %>%
  pivot_wider(names_from = amplification,values_from = type)

test22 <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  mutate(test2_1 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,1],
         test2_2 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,2],
         test2_3 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,3]) %>%
  filter(str_detect(test2_2,"amp")) %>%
  select(Participant,test2_2) %>%
  mutate(amplification = str_replace_all(test2_2,"_amp",""),type = "Amplication") %>% 
  select(Participant,amplification,type) %>%
  pivot_wider(names_from = amplification,values_from = type)

test23 <- df1 %>% select(Participant,ends_with('.pathway.alteration')) %>%
  mutate(test2_1 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,1],
         test2_2 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,2],
         test2_3 = str_split_fixed(PIK3CA.pathway.alteration,"\\|",3)[,3]) %>%
  filter(str_detect(test2_3,"amp")) %>%
  select(Participant,test2_3) %>%
  mutate(amplification = str_replace_all(test2_3,"_amp",""),type = "Amplication") %>%
  select(Participant,amplification,type) %>%
  pivot_wider(names_from = amplification,values_from = type)




amp <- test11 %>%full_join(test12) %>%
  full_join(test21) %>%
  full_join(test22) %>%
  full_join(test23)
amp

