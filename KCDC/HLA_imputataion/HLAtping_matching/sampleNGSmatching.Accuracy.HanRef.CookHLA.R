library(stringr)
library(tidyverse)

setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/cookHLA_han/acc/")
df <-read.table("JG.Han.HG18.HLAimp.MHC.HLA_IMPUTATION_OUT_withheader_511sample_fd.cmp_Nomencleaner.fd.txt",header = T)
head(df)

a <- df %>% summarise(across(colnames(df)[-1],sum))
for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
  a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
}
a
a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum))
a


df <-read.table("JG.Han.HG18.HLAimp.MHC.HLA_IMPUTATION_OUT_withheader_511sample_fd.cmp_Nomencleaner.fdvstd_td.txt",header = T)
head(df)

a <- df %>% summarise(across(colnames(df)[-1],sum))
for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
  a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
}
a
a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum))
a



library(stringr)
library(tidyverse)

setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/cookHLA_han/520/")
df <-read.table("KRKD.HLAimp_cookHLA.HANref.hg18.520sample.MHC.HLA_IMPUTATION_OUT_withheader.cmp_Nomencleaner.fd.txt",header = T)
head(df)

a <- df %>% summarise(across(colnames(df)[-1],sum))
for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
  a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
}
a
a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum))
a
b <- a

df <-read.table("KRKD.HLAimp_cookHLA.HANref.hg18.520sample.MHC.HLA_IMPUTATION_OUT_withheader.cmp_Nomencleaner.fdvstd.txt",header = T)
head(df)

a <- df %>% summarise(across(colnames(df)[-1],sum))
for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
  a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
}
a
a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum))
c <- rbind(a,b)



df <- read.table("KRKD.HLAimp_cookHLA.HANref.hg18.520sample.MHC.HLA_IMPUTATION_OUT_withheader.cmp_Nomencleaner.fd_missINFO.txt",header = T)
head(df)
wrong_list <- df %>% select(gene,Real) %>% unique()
ref <- read.table("../../../../HLA_type_frequency.txt",header = T)
head(ref)

ref %>% filter(digit == "4digit") %>% filter(HLA_type %in% wrong_list$Real) %>%
  mutate(HLA_gene = str_split_fixed(HLA_gene,"-",2)[,2]) %>% 
  mutate(freq = as.numeric(freq)) %>%
  ggplot(aes(x=HLA_gene,y=freq,fill=HLA_gene)) +
  geom_boxplot() + 
  theme(legend.position = 'none')
  
  
head(c)



df <-read.table("JG.Han.HG18.HLAimp.MHC.HLA_IMPUTATION_OUT_withheader_511sample_fd.cmp_Nomencleaner.fd.txt",header = T)
head(df)

#a <- df %>% summarise(across(colnames(df)[-1],sum))
a <- df
for (gene in c("HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1")) {
  a[,gene] <- a[,paste0(gene,".match")]/(a[,paste0(gene,".match")] + a[,paste0(gene,".wrong")])
}
a
a <- a %>% mutate("overall" = match_Sum/(match_Sum + wrong_Sum))
head(a)


a %>% select(c("ID","HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1","overall")) %>% #head()
  pivot_longer(cols = 2:10,names_to = "Gene",values_to = "Accuracy") %>% #count(Gene)
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>% #count(Gene)#head()
  ggplot(aes(x=Gene,y=Accuracy,fill=Gene)) +
  geom_boxplot()
  #geom_violin()
