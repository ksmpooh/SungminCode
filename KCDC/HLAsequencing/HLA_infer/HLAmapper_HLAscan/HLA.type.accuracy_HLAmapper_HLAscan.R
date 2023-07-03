#### HLA mapper and HLA scan 
#### HLA infor accuracy

library(stringr)
library(tidyverse)

setwd("~/Desktop/KCDC/HLA_seq/HLAinfer/HLA_scan/")


################
#flist = grep(list.files("./"),pattern = "missINFO.txt", invert=TRUE, value=TRUE)
flist = grep(list.files("./"),pattern = "cmp_Nomencleaner", invert=FALSE, value=TRUE)
flist = grep(flist,pattern = "missINFO", invert=TRUE, value=TRUE)
flist
out = NULL
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

out %>% mutate(filename = str_replace_all(filename,".txt","")) %>%# head()
  mutate('digit' = ifelse(grepl(pattern = "td",filename),"2","4")) %>%
  mutate("CV" = str_split_fixed(filename,"\\.",3)[,2]) %>% #head()
  mutate("Ref" = str_split_fixed(filename,"\\.",4)[,4]) -> out


fd <-read.table("merge.df.cmp_Nomencleaner.fd.txt",header = T)
head(fd)
fd %>% mutate("Accuracy" = match_Sum/(wrong_Sum+match_Sum)*100) %>% 
  ggplot(aes(y=Accuracy)) +
  geom_boxplot()
head(fd)

grep(".wrong",colnames(fd))
fd %>% select(grep(".wrong",colnames(fd))) %>%
  pivot_longer(1:8,values_to = "Count",names_to = 'Gene') %>% 
  mutate(Gene = str_replace_all(Gene,".wrong","")) %>%
  mutate(Gene = str_replace_all(Gene,"HLA_","")) %>%
  ggplot(aes(x=Gene,y=Count,color=Gene)) +
  geom_boxplot()
  
