## 20230125 HLA gene mapping coverage

library(tidyverse)
library(stringr)
setwd("~/Desktop/KCDC/HLA_seq/02.bam.stat_HLAgene/")

wDir <- "longread/coverage_hg19/"
a <- list.files(wDir)
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
df
a[1]
str_split_fixed(a[1],".coverage_",2)
str_split_fixed(a[1],"_HLA-",2)
str_split_fixed(str_split_fixed(a[1],"_HLA-",2)[,1],".coverage_",2)
#추가 해야할 것 : sampleID, build, contig, HLA gene, 
# sampleID
str_split_fixed(a[1],"\\.",7)
#       [,1]  [,2]       [,3]  [,4]  : ID  [,5]    reference file     [,6]  [,7]   contig, HLA-A                             
#[1,] "HLA" "Longread" "Seq" "NIH19KT0247" "hg19_HLAregion_mapped" "bam" "coverage_6_28477797_33448354_HLA-A"



## long hg19
wDir <- "longread/coverage_hg19/"
a <- list.files(wDir)
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> long_hg19

## long hg38
wDir <- "longread/coverage_hg38_HLAregion/"
a <- list.files(wDir)
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> long_hg38


## long hg38 with alt
wDir <- "longread/coverage_hg38_HLAregion_withALT/"
a <- list.files(wDir)
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> long_hg38_alt

head(long_hg19)
head(long_hg38)
head(long_hg38_alt)

long_hg19 %>% rbind(long_hg38) %>% rbind(long_hg38_alt) %>% #head()
  ggplot()


long_hg38_alt %>% select(-startpos,-endpos) %>% #head()
  ggplot(aes(y=meandepth,x=contig,color=contig)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 2)

head(long_hg38)
head(long_hg19)
long_hg38_alt %>% count(rname,contig)
long_hg38 %>% count(rname,contig)

long_hg38 %>% mutate(rname = "(only)chr6_28510020_33480577") %>%
  rbind(long_hg19 %>% mutate(rname = "(hg19)6_28477797_33448354")) %>%
  rbind(long_hg38_alt) %>% select(-startpos,-endpos) %>% #head()
  ggplot(aes(y=meandepth,x=rname,color=rname)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 2) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

long <- long_hg38 %>% mutate(rname = "(only)chr6_28510020_33480577") %>%
  rbind(long_hg19 %>% mutate(rname = "(hg19)6_28477797_33448354")) %>%
  rbind(long_hg38_alt) %>% mutate(HLAseq = "longread")
head(long)


######### short-read

## long hg19
wDir <- "shortread/coverage_hg19/"
a <- list.files(wDir)
#a <- a[grepl("dedup",a)]
a <- a[grepl("sorted.bam",a)]
a
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> short_hg19

## long hg38
wDir <- "shortread/coverage_hg38_HLAregion/"
a <- list.files(wDir)
#a <- a[grepl("dedup",a)]
a <- a[grepl("sorted.bam",a)]
a
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> short_hg38

head(short_hg38)


## long hg38 alt
wDir <- "shortread/coverage_hg38_HLAregion_withALT/"
a <- list.files(wDir)
#a <- a[grepl("dedup",a)]
a <- a[grepl("sorted.bam",a)]
a
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> short_hg38_alt

head(short_hg38_alt)


short <- short_hg38 %>% mutate(rname = "(only)chr6_28510020_33480577") %>%
  rbind(short_hg19 %>% mutate(rname = "(hg19)6_28477797_33448354")) %>%
  rbind(short_hg38_alt) %>% mutate(HLAseq = "shortread")
head(short)

short %>% ggplot(aes(y=meandepth,x=rname,fill=rname)) +
  geom_boxplot() +
  #facet_wrap(~HLAgene,nrow = 2) +
  facet_wrap(~HLAgene,nrow = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#### dedup
## long hg19
wDir <- "shortread/coverage_hg19/"
a <- list.files(wDir)
a <- a[grepl("dedup",a)]
#a <- a[grepl("sorted.bam",a)]
a
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> short_hg19_dedup
head(short_hg19_dedup)

## long hg38
wDir <- "shortread/coverage_hg38_HLAregion/"
a <- list.files(wDir)
a <- a[grepl("dedup",a)]
#a <- a[grepl("sorted.bam",a)]
a
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> short_hg38_dedup

head(short_hg38_dedup)


## short hg38 alt
wDir <- "shortread/coverage_hg38_HLAregion_withALT/"
a <- list.files(wDir)
a <- a[grepl("dedup",a)]
#a <- a[grepl("sorted.bam",a)]
a
df <- read.table(paste0(wDir,a[1]))
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
colnames(df) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
df$filename <- a[1] 
for (i in a[2:length(a)]) {
  tmp <- read.table(paste0(wDir,i))
  colnames(tmp) <- c('rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq')
  tmp$filename <- i
  #print(i)
  df <- rbind(df,tmp)
}
head(df)

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> short_hg38_alt_dedup

head(short_hg38_alt)


short_dedup <- short_hg38_dedup %>% mutate(rname = "(only)chr6_28510020_33480577") %>%
  rbind(short_hg19_dedup %>% mutate(rname = "(hg19)6_28477797_33448354")) %>%
  rbind(short_hg38_alt_dedup) %>% mutate(HLAseq = "shortread")
head(short)

short_dedup %>% ggplot(aes(y=meandepth,x=rname,fill=rname)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 2) +
  #facet_wrap(~HLAgene,nrow = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




### plot

long %>% ggplot(aes(y=meandepth,x=rname,fill=rname)) +
  geom_boxplot() +
  #facet_wrap(~HLAgene,nrow = 2) +
  facet_wrap(~HLAgene,nrow = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

short_dedup %>% ggplot(aes(y=meandepth,x=rname,fill=rname)) +
  geom_boxplot() +
  #facet_wrap(~HLAgene,nrow = 2) +
  facet_wrap(~HLAgene,nrow = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

short_dedup %>% ggplot(aes(y=meandepth,x=rname,fill=rname)) +
  geom_boxplot() +
  #facet_wrap(~HLAgene,nrow = 2) +
  facet_wrap(~HLAgene,nrow = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")
