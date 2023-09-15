#### HLA type
library(tidyverse)
library(stringr)
library(readxl)

setwd("~/Desktop/KCDC/????????????????????????????????????/??????????????????????????????/")
ref <- read_xlsx("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/HLA.typing.Final.result_529sample.xlsx")
ref <- ref %>% select(Sample,ID)
head(ref)
ind <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Final_520sample.index.txt",header = T)
head(ind)
a <- read_xlsx("2019/all_sample_genotypes.xlsx",sheet = 3,skip = 1)
head(a);colnames(a)

#a %>% select(1,grep("Allele",colnames(a))) -> t
a %>% select(1,grep("Allele",colnames(a))) %>% 
  merge(ref,by.x='...1',by.y="Sample") %>% filter(ID %in% ind$KBAID) %>% select(-...1)-> a

colnames(a) <- c("A_1","A_2","B_1","B_2","C_1","C_2","DRB1_1","DRB1_2","DRB3_1","DRB3_2","DQA1_1","DQA1_2","DQB1_1","DQB1_2","DPA1_1","DPA1_2","DPB1_1","DPB1_2","ID")

b <- read_excel("2020_NGgene???????????????/3.HLA genotyping_2020_265 samples 11 locus_????????????????????????.xlsx",sheet = 5)


b %>% merge(ref) %>% filter(ID %in% ind$KBAID) %>% select(-Sample) -> b
colnames(b) <- c("A_1","A_2","B_1","B_2","C_1","C_2","DRB1_1","DRB1_2","DRB3_1","DRB3_2","DRB4_1","DRB4_2","DRB5_1","DRB5_2","DQA1_1","DQA1_2","DQB1_1","DQB1_2","DPA1_1","DPA1_2","DPB1_1","DPB1_2","ID")
head(b)

a %>% select(ID,grep("_1",colnames(a))) %>% pivot_longer(2:10) %>% 
  mutate(gene = str_split_fixed(name,"_",2)[,1]) %>%
  mutate(value1 = str_split_fixed(value,"/",2)[,1]) %>% #head()
  select(-name,-value) -> a.1

a %>% select(ID,grep("_2",colnames(a))) %>% pivot_longer(2:10) %>% 
  mutate(gene = str_split_fixed(name,"_",2)[,1]) %>%
  mutate(value2 = str_split_fixed(value,"/",2)[,1]) %>% #head()
  select(-name,-value) -> a.2

b %>% select(ID,grep("_1",colnames(b))) %>% pivot_longer(2:12) %>% 
  mutate(gene = str_split_fixed(name,"_",2)[,1]) %>%
  mutate(value1 = str_split_fixed(value,"/",2)[,1]) %>% #head()
  select(-name,-value) -> b.1

b %>% select(ID,grep("_2",colnames(b))) %>% pivot_longer(2:12) %>% 
  mutate(gene = str_split_fixed(name,"_",2)[,1]) %>%
  mutate(value2 = str_split_fixed(value,"/",2)[,1]) %>% #head()
  select(-name,-value) -> b.2

head(a.1);dim(a.1)
head(a.2);dim(a.2)
head(b.1);dim(b.1)
head(b.2);dim(b.2)

a.1 %>% left_join(a.2) -> a.0
b.1 %>% left_join(b.2) -> b.0
dim(a.0)
dim(b.0)

head(b.0) 
b.0 %>% filter(!gene %in% c("DRB3","DRB4","DRB5")) %>% filter(value1==value2) 

a.0 %>% left_join(b.0) %>% count(gene)
  filter(value1 ==".",value2 == ".") %>% filter(gene %in% c("DPA1","DRB1"))

a.0 %>% rbind(b.0) %>% filter(!gene %in% c("DRB3","DRB4","DRB5")) %>% #head()
  mutate(value1 = ifelse(value1 == ".",value2,value1),value2 = ifelse(value2 == ".",value1,value2)) %>%  #dim()
  pivot_longer(3:4) -> g_type

g_type %>% count(gene)
g_type %>% group_by(gene) %>% count(value) %>% filter(!value %in% c(".","None")) %>%
  mutate(Frequency = prop.table(n)) %>% select(-n) %>% writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/99.forPaper/HLA.type.520samples.Ggroup.frequency.xlsx")

g_type %>% group_by(gene) %>% count(value) %>% filter(!value %in% c(".","None")) %>%
  mutate(Frequency = prop.table(n)) -> g_freq
g_freq %>% group_by(gene) %>% top_n(1,n) %>% select(gene,value) %>% rename(max = value)-> g_max_value

head(g_max_value)


a.0 %>% rbind(b.0) %>% filter(!gene %in% c("DRB3","DRB4","DRB5")) %>% #count(gene)
  mutate(value1 = ifelse(value1 == ".",value2,value1),value2 = ifelse(value2 == ".",value1,value2)) %>% #head()
  mutate(value1 = ifelse(value1 == "None",value2,value1),value2 = ifelse(value2 == "None",value1,value2)) %>% #head()
  pivot_longer(3:4) %>% left_join(g_max_value) %>%
  group_by(gene) %>% 
  mutate(value = ifelse(value %in% c("None","."),max,value)) %>% #filter(ID %in% c("NIH19KT2541") & gene %in% c("A"))#filter(ID %in% c("NIH19KT0003","NIH19KT5872","NIH19KT6013","NIH19KT6050","NIH19KT2541") & gene %in% c("A","DRB1","DPA1","DBP1"))
  mutate(value = str_split_fixed(value,"\\*",2)[,2]) %>%# mutate(value = str_replace_all(value,"G",'g')) %>%
  mutate(name = paste0(gene,"_",name)) %>% ungroup() %>%select(-gene,-max) %>% #head()
  pivot_wider(names_from = "name",values_from = "value") -> g_forhalp

g_forhalp %>% dim()

head(g_forhalp)

colnames(g_forhalp) <- c("id","A","A","B","B","C","C","DRB1","DRB1","DQA1","DQA1","DQB1","DQB1","DPA1","DPA1","DPB1","DPB1")

write.table(g_forhalp,"~/Desktop/KCDC/HLAimputation/99.forPaper/HLA.type.520samples.Ggroup.forHaplot_O_mat.txt",col.names = T,row.names = F,sep = "\t",quote = F)




head(a.)