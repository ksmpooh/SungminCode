### Q30 bam stat compare

## hg19_whole vs hg19_HLAregion_target
## hg38_whole vs hg38_HLAregion_target + alt

library(tidyverse)
library(stringr)
library(grid)
library(cowplot)

setwd("~/Desktop/KCDC/HLA_seq/02.bam.stat_HLAgene/")

wDir <- "longread/Q30/coverage_hg19_gene/"
a <- list.files(wDir)
df <- NULL
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
for (i in a[1:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> long_hg19_HLAtarget_Q30



wDir <- "longread/coverage_hg19_Q20/"
a <- list.files(wDir)
df <- NULL
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
for (i in a[1:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> long_hg19_HLAtarget_Q20

head(long_hg19_HLAtarget_Q30)
head(long_hg19_HLAtarget_Q20)


long_hg19_HLAtarget_Q30$Q <- "Q30"
long_hg19_HLAtarget_Q20$Q <- "Q20"

long_hg19_HLAtarget_Q30 %>% rbind(long_hg19_HLAtarget_Q20) %>% 
  select(ID,Q,HLAgene,meandepth,coverage) %>% #head()
  ggplot(aes(x=Q,y=meandepth,fill=Q)) + 
  geom_boxplot() + 
  facet_grid(~HLAgene)


#### hg38

wDir <- "longread/Q30/coverage_hg38_alt_gene/"
a <- list.files(wDir)
df <- NULL
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
for (i in a[1:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) ->  long_hg38_HLAtarget_Q30



wDir <- "longread/Q30/coverage_hg38_whole_gene/"
a <- list.files(wDir)
df <- NULL
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
for (i in a[1:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> long_hg38_whole_Q30

head(long_hg38_HLAtarget_Q30)
head(long_hg38_whole_Q30)


long_hg38_HLAtarget_Q30 %>% rbind(long_hg38_whole_Q30) %>% count(contig)
  ggplot(aes(x=type,y=meandepth,fill=type)) + 
  geom_boxplot() +
  facet_grid(~HLAgene)





