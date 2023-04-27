library(readxl)
library(tidyverse)
#df <-read.csv("~/Downloads/KORV1_1_na35_annot_open.csv")
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/")


'mutate(gene = str_replace_all(gene,"NGS_","")) %>%
  mutate(gene = str_replace_all(gene,"\\.1","")) %>%
  mutate(gene = str_replace_all(gene,"\\.2","")) %>% #head()
  '

df_td <-read.table("test/michigan/02.processing/KBA.520sample.michiganHLAimp_td.txt",header = T)
df_fd <-read.table("test/michigan/02.processing/KBA.520sample.michiganHLAimp_fd.txt",header = T)
head(df_td)
head(df_fd)

df_td %>% pivot_longer(2:17,names_to = 'gene',values_to = 'imp_td_type') %>% head()
df_td %>% pivot_longer(grep("\\.1",colnames(df_td)),names_to = 'gene',values_to = 'imp_td_type1') %>%
  mutate(gene = str_replace_all(gene,"\\.1","")) %>%
  select(ID,gene,imp_td_type1) -> df_td_1
df_td %>% pivot_longer(grep("\\.2",colnames(df_td)),names_to = 'gene',values_to = 'imp_td_type2') %>% 
  mutate(gene = str_replace_all(gene,"\\.2","")) %>%
  select(ID,gene,imp_td_type2) -> df_td_2
  #pivot_longer(colnames(df_td)[grepl("\\.2",colnames(df_td))],names_to = 'gene2',values_to = 'imp_td_type2') %>% head()
head(df_td_1)
head(df_td_2)

df_td_1 %>% full_join(df_td_2) -> df_td
  
df_fd %>% pivot_longer(2:17,names_to = 'gene',values_to = 'imp_td_type') %>% head()
df_fd %>% pivot_longer(grep("\\.1",colnames(df_fd)),names_to = 'gene',values_to = 'imp_fdtd_type1') %>%
  mutate(gene = str_replace_all(gene,"\\.1","")) %>%
  mutate(imp_fdtd_type1 = str_split_fixed(imp_fdtd_type1,":",2)[,1]) %>%
  select(ID,gene,imp_fdtd_type1) -> df_fd_1

df_fd %>% pivot_longer(grep("\\.2",colnames(df_fd)),names_to = 'gene',values_to = 'imp_fdtd_type2') %>% 
  mutate(gene = str_replace_all(gene,"\\.2","")) %>%
  mutate(imp_fdtd_type2 = str_split_fixed(imp_fdtd_type2,":",2)[,1]) %>%
  select(ID,gene,imp_fdtd_type2) -> df_fd_2

head(df_fd_1)
head(df_fd_2)

df_fd_1 %>% full_join(df_fd_2) -> df_fd

head(df_fd)

df <- df_td %>% full_join(df_fd)
head(df)

#  select(-gene1)


ref <-read.table("HLA.type.result.8genes.merged.4digit_529sample_forMAKEreference.txt",header = T)
head(ref)

ref %>% select(-FID,-pID,-mID,-SEX,-PHENO) %>% 
  rename("ID" = IID) -> ref0

head(ref0)
  
ref0 %>% pivot_longer(grep("\\.1",colnames(ref0)),names_to = 'gene',values_to = 'NGS_type1') %>%
  mutate(gene = str_replace_all(gene,"\\.1","")) %>%
  mutate(NGS_type1 = str_split_fixed(NGS_type1,":",2)[,1]) %>%
  select(ID,gene,NGS_type1) -> ref1

ref0 %>% pivot_longer(grep("\\.2",colnames(ref0)),names_to = 'gene',values_to = 'NGS_type2') %>% 
  mutate(gene = str_replace_all(gene,"\\.2","")) %>%
  mutate(NGS_type2 = str_split_fixed(NGS_type2,":",2)[,1]) %>%
  select(ID,gene,NGS_type2) -> ref2
head(ref1)
head(ref2)

ref1 %>% full_join(ref2) %>% mutate(gene = str_replace_all(gene,"NGS","HLA")) %>% 
  filter(ID %in% df$ID)-> ref

head(df)
head(ref)

m_ref1 <- read.table("test/michigan/03.allele.matching/KBA.520sample.michiganHLAimp_fd.cmp_Nomencleaner.fdvstd.txt",header = T)
m_ref2 <- read.table("test/michigan/03.allele.matching/KBA.520sample.michiganHLAimp_td.cmp_Nomencleaner.txt",header = T)
head(m_ref)
colnames(m_ref)[grep('wrong',colnames(m_ref))]
head(m_ref1)
m_ref1 %>% select(ID,grep('wrong',colnames(m_ref1))) %>%
  pivot_longer(colnames(m_ref1)[grep('wrong',colnames(m_ref1))],names_to = 'gene',values_to = 'wrong_fdtd') %>% 
  mutate(gene=str_replace_all(gene,".wrong","")) %>% filter(gene != "wrong_Sum") -> m_ref1


m_ref2 %>% select(ID,grep('wrong',colnames(m_ref2))) %>%
  pivot_longer(colnames(m_ref2)[grep('wrong',colnames(m_ref2))],names_to = 'gene',values_to = 'wrong_td') %>% 
  mutate(gene=str_replace_all(gene,".wrong","")) %>% filter(gene != "wrong_Sum") -> m_ref2
  #head()

head(m_ref1)
head(m_ref2)

df %>% full_join(ref) %>% right_join(m_ref1 %>% full_join(m_ref2)) -> out# -> writexl::write_xlsx("michigan.tdfd.chekc.xlsx")
writexl::write_xlsx(out,"michigan.tdfd.chekc.xlsx")
