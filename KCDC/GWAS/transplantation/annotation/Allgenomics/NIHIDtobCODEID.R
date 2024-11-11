library(tidyverse)

setwd("~/Desktop/KCDC/transplantation/allogenomic/gene/")


a <- read_table("KR.KD.Main_Category.alleleMatching01.Score_Sum.txt")
b <- read_table("KR.KD.MHC_I_each_gene.alleleMatching01.Score_Sum.txt")
c <- read_table("KR.KD.MHC_II_each_gene.alleleMatching01.Score_Sum.txt")
d <- read_table("KR.KD.transmembrane_each_gene.alleleMatching01.Score_Sum.txt")

ref <- read_table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt")
head(ref)
head(a)
head(b)
head(c)
head(d)



ref %>% rename(KBA_ID.KD = KBA_ID) %>% left_join(a)
a %>% left_join(ref %>% rename(KBA_ID.KR = KBA_ID,bCODE_R = bCODE)) %>%
  left_join(ref %>% rename(KBA_ID.KD = KBA_ID,bCODE_D = bCODE)) %>% select(ncol(a)+1,ncol(a)+2,3:ncol(a)) -> a
head(a)

colnames(a)[3:ncol(a)] <- str_replace(colnames(a)[3:ncol(a)],"_Score_SUM","")
a %>% select(1,2,3,6,7,4,5,8) -> a


head(d)
d[str_detect(colnames(d),"HLA")] -> d_outlist
d %>% select(-colnames(d_outlist)) %>% select(DOB_Score_SUM)

'''
2 HLA-DOA
2 HLA-DOB
2 HLA-DPA1
2 HLA-DQA1
2 HLA-DQA2
2 HLA-DQB1
2 HLA-DQB2
2 HLA-DRA
2 HLA-DRB1
'''


head(d)


d %>% left_join(ref %>% rename(KBA_ID.KR = KBA_ID,bCODE_R = bCODE)) %>%
  left_join(ref %>% rename(KBA_ID.KD = KBA_ID,bCODE_D = bCODE)) %>% select(ncol(a)+1,ncol(a)+2,3:ncol(a)) -> a
head(a)

colnames(a)[3:ncol(a)] <- str_replace(colnames(a)[3:ncol(a)],"_Score_SUM","")
head(b)

b %>% mutate(MHC_I = rowSums(across(3:last_col()))) %>% select(KBA_ID.KD,KBA_ID.KR,MHC_I) %>% left_join(a) %>% count(MHC_I == MHC_I_Score_SUM)
c %>% 
  #mutate(MHC_I = rowSums(across(3:last_col()))) %>% 
  select(KBA_ID.KD,KBA_ID.KR,MHC_II_Score_SUM) %>% rename(MHC_II = MHC_II_Score_SUM) %>% #head()
  left_join(a) %>% count(MHC_II == MHC_II_Score_SUM)

head(c)


head(a)
head(c)
head(d)
head(a)
d[,c(1,2)]





d %>%  select(3000:3004)
c %>% left_join(d) %>% select(3000:3004)
b %>% left_join(c) -> o0
d %>% left_join(o0) %>% select(3540:3545)
o0 %>% left_join(d) 

o0 %>% select(KBA_ID.KD,KBA_ID.KD) %>% m
o0 %>% count(KBA_ID.KD %in% d$KBA_ID.KD)

b %>% left_join(c) %>% #head()
  full_join(d %>% select(-colnames(d_outlist))) %>% #select(3000:3004)
  left_join(ref %>% rename(KBA_ID.KR = KBA_ID,bCODE_R = bCODE)) %>%
  left_join(ref %>% rename(KBA_ID.KD = KBA_ID,bCODE_D = bCODE)) -> o1

head(o1)
head(a)
#head(a)
o1 %>% select(ncol(o1)-1,ncol(o1),3:(ncol(o1)-2)) -> o1
colnames(o1)[3:ncol(o1)] <- str_replace(colnames(o1)[3:ncol(o1)],"_Score_SUM","")
head(o1)
a %>% left_join(o1) %>% select(-transmembrane) %>%
  write.table("AMS.score.byGene.MHC_transmemebrane.txt",col.names = T,row.names = F,quote = F,sep = "\t")



