## 20230125 HLA gene mapping coverage

library(tidyverse)
library(stringr)
library(grid)
library(cowplot)
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

## short hg19
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

## short hg38
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
head(short_dedup)

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


### by exon 20230213
library(tidyverse)
library(stringr)
setwd("~/Desktop/KCDC/HLA_seq/02.bam.stat_HLAgene/")

#long 19 exon
wDir <- "longread/coverage_hg19_exon/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> long_hg19_exon

## long hg38
wDir <- "longread/coverage_hg38_HLAregion_exon/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5],'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> long_hg38_exon
#head(long_hg38_exon)



## short exon

wDir <- "shortread/coverage_hg19_exon/"
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


df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> short_hg19_exon
head(short_hg19_exon)

## short hg38 exon
wDir <- "shortread/coverage_hg38_HLAregion_exon/"
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

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6],'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> short_hg38_exon


head(short_hg19_exon)
head(short_hg38_exon)

head(long_hg19_exon)
head(long_hg38_exon)

## plot


long_hg19_exon %>% ggplot(aes(y=meandepth,x=exon,fill=exon)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 2) +
#  facet_wrap(~HLAgene,nrow = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

long_hg38_exon %>% ggplot(aes(y=meandepth,x=exon,fill=exon)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 2) +
  #  facet_wrap(~HLAgene,nrow = 1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

long_hg19_exon %>% rbind(long_hg38_exon) %>% #head()
  ggplot(aes(y=meandepth,x=exon,fill=type)) +
  geom_boxplot() + 
  facet_wrap(~HLAgene,nrow = 2) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")


long_hg19_exon %>% rbind(long_hg38_exon) %>% head()
short_hg19_exon %>% rbind(short_hg38_exon) %>% #head()
  ggplot(aes(y=meandepth,x=exon,fill=type)) +
  geom_boxplot() + 
  facet_wrap(~HLAgene,nrow = 2) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")

'
seq       build meandepth meanbaseq
<chr>     <chr>     <dbl>     <dbl>
  1 Longread  hg19       132.      87.4
2 Longread  hg38       132.      87.4
3 Shortread hg19       398.      36.4
4 Shortread hg38       398.      36.4
'

## 5M


library(tidyverse)
library(stringr)
setwd("~/Desktop/KCDC/HLA_seq/02.bam.stat_HLAgene/")


wDir <- "longread/coverage_hg38_whole/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5]) %>%
  select(-filename) -> long_hg38_whole

head(long_hg38_whole)

wDir <- "longread/coverage_hg19_whole/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,5]) %>%
  select(-filename) -> long_hg19_whole

## short hg19
wDir <- "shortread/coverage_hg19_whole/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6]) %>% 
  select(-filename) -> short_hg19_whole

wDir <- "shortread/coverage_hg38_whole/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,4],'type'=str_split_fixed(filename,"\\.",7)[,6]) %>% 
  select(-filename) -> short_hg38_whole


#GRch38
head(long_hg19_whole)
head(long_hg38_whole)
head(short_hg38_whole)
head(short_hg19_whole)


long_hg19_whole %>% rbind(long_hg38_whole) %>% mutate("seq" = "Longread") %>%
  rbind(short_hg19_whole %>% rbind(short_hg38_whole) %>% mutate("seq" = "Shortread")) %>%
  mutate('build' = ifelse(str_detect(type,"hg19"),"hg19","hg38")) -> whole_seq

whole_seq %>% filter(seq == "Longread") %>% 
  ggplot(aes(y=meandepth,x=build,fill=build)) +
  geom_boxplot()

head(whole_seq)
whole_seq %>% group_by(seq,build) %>% #head()
  summarise(meandepth = mean(meandepth),
            coverage = mean (coverage),
            meanmapq = mean (meanmapq),
            numreads = mean (numreads))


##########mapping depth with coverage
title <- function(x){
  ggdraw(x) + 
    draw_label(
      x,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 7)
    )
}


head(long_hg19_exon)
long_hg19_exon %>% ggplot() +
  geom_boxplot(aes(y=meandepth,x=exon,fill=exon)) +
  geom_hline(yintercept= 132, linetype='dashed', color='red') +
#  geom_line(aes(y=mean(coverage),x=exon)) + 
  #facet_wrap(~HLAgene,nrow = 2) +
  facet_wrap(~HLAgene,nrow = 1) +
  labs(title = element_text("Longread")) + 
  theme(#axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        legend.position = "none") -> p1
p1

head(short_hg19_exon)
short_hg19_exon %>% ggplot() +
  geom_boxplot(aes(y=meandepth,x=exon,fill=exon)) +
  geom_hline(yintercept= 398, linetype='dashed', color='red') +
  #  geom_line(aes(y=mean(coverage),x=exon)) + 
  #facet_wrap(~HLAgene,nrow = 2) +
  facet_wrap(~HLAgene,nrow = 1) +
  labs(title = element_text("Shortread")) + 
  theme(#axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    legend.position = "none") -> p2
p2
plot_grid(title("Mean Depth by Exon: hg19"),
  p1, p2,
  labels = c("","A","B"),
  rel_heights = c(0.1, 1,1), ncol = 1)




head(long_hg19)
long_hg19 %>% rbind(long_hg38) %>% ggplot() +
  geom_boxplot(aes(y=meandepth,x=HLAgene,fill=HLAgene)) +
  facet_wrap(~type,ncol = 2) +
  theme(#axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") -> p1

long_hg19 %>% rbind(long_hg38) %>% ggplot() +
  geom_boxplot(aes(y=coverage,x=HLAgene,fill=HLAgene)) +
  facet_wrap(~type,ncol = 2) +
  
    theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p2

long_hg19 %>% rbind(long_hg38) %>% ggplot() +
  geom_boxplot(aes(y=meanmapq,x=HLAgene,fill=HLAgene)) +
  facet_wrap(~type,ncol = 2) +
  
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p3

long_hg19 %>% rbind(long_hg38) %>% ggplot() +
  geom_boxplot(aes(y=numreads,x=HLAgene,fill=HLAgene)) +
  facet_wrap(~type,ncol = 2) +
  
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p4

plot_grid(
  p1, p2, p3,p4,
  labels = "AUTO", ncol = 1)



head(long_hg19)
head(short_hg19_dedup)

long_hg19 %>% rbind(short_hg19_dedup) %>% #head()
  mutate("Seq" = ifelse(str_detect(type,'align'),"Shortread","Longread")) %>%
  mutate("averagedepth" = ifelse(Seq=="Shortread","398","132"))-> df

long_hg19 %>% ggplot() +
  geom_boxplot(aes(y=meandepth,x=HLAgene,fill=HLAgene)) +
  geom_hline(yintercept= 132, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p1

long_hg19 %>% ggplot() +
  geom_boxplot(aes(y=coverage,x=HLAgene,fill=HLAgene)) +
#  geom_hline(yintercept= 132, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p2

long_hg19 %>% ggplot() +
  geom_boxplot(aes(y=meanmapq,x=HLAgene,fill=HLAgene)) +
  #geom_hline(yintercept= 132, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p3

long_hg19 %>% ggplot() +
  geom_boxplot(aes(y=numreads,x=HLAgene,fill=HLAgene)) +
  #geom_hline(yintercept= 132, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p4


title <- function(x){
  ggdraw(x) + 
    draw_label(
      x,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 7)
    )
}


plot_grid(p1, p2, p3,p4,
          labels = c("A","B","C","D"), nrow = 1) -> a1

plot_grid(title("Longread"),a1,
          ncol = 1,
          rel_heights = c(0.1, 2)
          ) ->a1


# 
# a1<- plot_grid(title("Longread"),
#   p1, p2, p3,p4,
#   labels = c("","A","B","C","D"), ncol = 1,
#   rel_heights = c(0.1, 1,1,1,1))
# a1


short_hg19_dedup %>% ggplot() +
  geom_boxplot(aes(y=meandepth,x=HLAgene,fill=HLAgene)) +
  geom_hline(yintercept= 398, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p1

short_hg19_dedup %>% ggplot() +
  geom_boxplot(aes(y=coverage,x=HLAgene,fill=HLAgene)) +
  #geom_hline(yintercept= 398, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p2

short_hg19_dedup %>% ggplot() +
  geom_boxplot(aes(y=meanmapq,x=HLAgene,fill=HLAgene)) +
  #geom_hline(yintercept= 398, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p3

short_hg19_dedup %>% ggplot() +
  geom_boxplot(aes(y=numreads,x=HLAgene,fill=HLAgene)) +
  #geom_hline(yintercept= 398, linetype='dashed', color='red') +
  theme(#axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p4


plot_grid(p1, p2, p3,p4,
          #labels = c("A","B","C","D"), 
          nrow = 1) -> a2

plot_grid(title("Shortread"),a2,
          ncol = 1,
          rel_heights = c(0.1, 2)
) ->a2


plot_grid(
  a1,a2,nrow=2
)

