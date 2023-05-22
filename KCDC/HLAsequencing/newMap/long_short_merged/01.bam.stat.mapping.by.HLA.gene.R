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
## short hg19
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


###### 20230407 merge mapping theme1,2 , short+long
library(tidyverse)
library(stringr)
library(grid)
library(cowplot)
setwd("~/Desktop/KCDC/HLA_seq/merge_align/bam.stat/")
#setwd("~/Desktop/KCDC/HLA_seq/merge_align/bam.stat/02_theme/")

wDir <- "01_theme/coverage_hg19_gene/"
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



## merged gene theme1
wDir <- "01_theme/coverage_hg19_gene/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,2],'type'="Theme 1",'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> theme1_hg19_gene
head(theme1_hg19_gene)


## merged gene theme2
wDir <- "02_theme/coverage_hg19_gene/"
#wDir <- "02_theme_dedup/coverage_hg19_gene/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 2",'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> theme2_hg19_gene
head(theme2_hg19_gene)

wDir <- "02_theme_dedup/coverage_hg19_gene/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 2 (dedup)",'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> theme2_dedup_hg19_gene
head(theme2_dedup_hg19_gene)

### 03 theme gene
wDir <- "03_theme/coverage_hg19_gene/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 3",'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> theme3_hg19_gene
head(theme3_hg19_gene)

wDir <- "03_theme_dedup/coverage_hg19_gene/"
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 3 (dedup)",'HLAgene'=str_split_fixed(filename,"_HLA-",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% select(-filename) -> theme3_dedup_hg19_gene
head(theme3_dedup_hg19_gene)


### 01 theme exon
wDir <- "01_theme/coverage_hg19_exon/"
a <- list.files(wDir)
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

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,2],'type'="Theme 1",'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> theme1_hg19_exon
head(theme1_hg19_exon)

### 02 theme exon
wDir <- "02_theme/coverage_hg19_exon/"
#wDir <- "02_theme_dedup/coverage_hg19_exon/"
a <- list.files(wDir)
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

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 2",'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> theme2_hg19_exon
head(theme2_hg19_exon)

wDir <- "02_theme_dedup/coverage_hg19_exon/"
a <- list.files(wDir)
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

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 2 (dedup)",'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> theme2_dedup_hg19_exon
head(theme2_dedup_hg19_exon)

#### theme 3 exon
wDir <- "03_theme/coverage_hg19_exon/"
#wDir <- "02_theme_dedup/coverage_hg19_exon/"
a <- list.files(wDir)
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

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 3",'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> theme3_hg19_exon
head(theme3_hg19_exon)

wDir <- "03_theme_dedup/coverage_hg19_exon/"
a <- list.files(wDir)
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

df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'="Theme 3 (dedup)",'HLAgene'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,1]) %>%
  mutate('exon'=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,2],"_exon",2)[,2]) %>%
  mutate("contig"=str_split_fixed(str_split_fixed(filename,"_HLA-",2)[,1],".coverage_",2)[,2]) %>% 
  select(-filename) -> theme3_dedup_hg19_exon
head(theme3_dedup_hg19_exon)



### HLA 5M whole
wDir <- "01_theme/coverage_hg19_HLAregion/"
a <- list.files(wDir)
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,2],'type'='Theme 1') %>% 
  select(-filename) -> theme1_hg19

head(theme1_hg19)


wDir <- "02_theme/coverage_hg19_HLAregion/"
#wDir <- "02_theme_dedup/coverage_hg19_HLAregion/"
a <- list.files(wDir)
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'='Theme 2') %>% 
  select(-filename) -> theme2_hg19

head(theme2_hg19)

wDir <- "02_theme_dedup/coverage_hg19_HLAregion/"
a <- list.files(wDir)
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'='Theme 2 (dedup)') %>% 
  select(-filename) -> theme2_dedup_hg19

head(theme2_dedup_hg19)

## theme 3 whole
wDir <- "03_theme/coverage_hg19_HLAregion/"
#wDir <- "02_theme_dedup/coverage_hg19_HLAregion/"
a <- list.files(wDir)
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'='Theme 3') %>% 
  select(-filename) -> theme3_hg19

head(theme3_hg19)

wDir <- "03_theme_dedup/coverage_hg19_HLAregion/"
a <- list.files(wDir)
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
df %>% mutate('ID'=str_split_fixed(filename,"\\.",7)[,3],'type'='Theme 3 (dedup)') %>% 
  select(-filename) -> theme3_dedup_hg19

head(theme3_dedup_hg19)




## plot

head(long_hg19_whole)
head(long_hg19)
head(long_hg19_exon)
head(short_hg19_whole)
head(short_hg19_exon)
head(short_hg19_dedup)


head(theme1_hg19)
head(theme1_hg19_exon)
head(theme1_hg19_gene)
head(theme2_hg19)
head(theme2_hg19_exon)
head(theme2_hg19_gene)
head(theme2_dedup_hg19)
head(theme2_dedup_hg19_gene)
head(theme2_dedup_hg19_exon)

head(theme3_hg19)
head(theme3_hg19_exon)
head(theme3_hg19_gene)
head(theme3_dedup_hg19)
head(theme3_dedup_hg19_gene)
head(theme3_dedup_hg19_exon)



## 5 M 
theme1_hg19 %>% rbind(theme2_hg19) -> merge_hg19
long_hg19_whole$type <- "Longread"
short_hg19_whole$type <- "Shortread"


head(long_hg19_whole)
colnames(long_hg19_whole)

#long_short_sum <- merge(long_hg19_whole,short_hg19_whole,all = T,by=ID)
head(long_short_sum)
long_hg19_whole %>% mutate(meandepth = meandepth + short_hg19_whole$meandepth) %>% 
  mutate(numreads = numreads + short_hg19_whole$numreads) %>% 
  #mutate(covbases = (covbases,short_hg19_whole$covbases)) %>% 
  mutate(coverage = (coverage + short_hg19_whole$coverage)/2) %>% 
  mutate(meanbaseq = (meanbaseq + short_hg19_whole$meanbaseq)/2) %>% 
  mutate(meanmapq = (meanmapq + short_hg19_whole$meanmapq)/2) %>%
  mutate(type = "Sum(Long+Short)") -> long_short_sum
  

long_hg19 %>% mutate(meandepth = meandepth + short_hg19_whole$meandepth) %>% 
  mutate(numreads = numreads + short_hg19_whole$numreads) %>% 
  mutate(type = "Sum(Long+Short)") -> long_short_sum_gene

head(long_short_sum_gene)







head(long_short_sum)
head(long_hg19_whole)
head(short_hg19_whole)
#head(short_hg19_)

theme1_hg19 %>% rbind(theme2_hg19) %>%
  rbind(long_short_sum) %>% #-> merge_hg19
  rbind(long_hg19_whole) %>%
  rbind(short_hg19_whole) -> merge_hg19

head(long_hg19_whole)

merge_hg19 %>%
  ggplot(aes(y=meandepth,x=type,fill=type)) +
  geom_boxplot() + 
  theme(axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none")-> p1
p1
  
merge_hg19 %>%
  ggplot(aes(y=numreads,x=type,fill=type)) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.title.x = element_blank(),
  legend.position = "none")-> p2


merge_hg19 %>%
  ggplot(aes(y=coverage,x=type,fill=type)) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p3

merge_hg19 %>%
  ggplot(aes(y=meanbaseq,x=type,fill=type)) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p4

merge_hg19 %>%
  ggplot(aes(y=meanmapq,x=type,fill=type)) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none") -> p5

merge_hg19 %>%
  ggplot(aes(y=meandepth,x=type,fill=type)) +
  geom_boxplot() +
  theme(legend.position = "bottom",
        legend.title = element_blank())-> p6

p6 <- get_legend(p6)
plot(p6)

plot_grid(p1, p2, p3,p4,p5,
          #labels = c("A","B","C","D"), 
          nrow = 1) -> a1


plot_grid(
  a1,p6,nrow=2,
  rel_heights = c(4, 0.5)
)


theme1_hg19 %>% rbind(theme2_hg19) %>%
  group_by(type) %>%
  summarise(meandepth = mean(meandepth),
            coverage = mean (coverage),
            meanmapq = mean (meanmapq),
            numreads = mean (numreads))

'
type    meandepth coverage meanmapq  numreads
<chr>       <dbl>    <dbl>    <dbl>     <dbl>
1 Theme 1      529.     99.9     49.4 22802959.
2 Theme 2      603.     99.9     43.9 26853313.
'



### 유전자별 결과 비교
head(long_hg19_whole)
head(long_hg19)

head(short_hg19_whole)
head(short_hg19_dedup)


head(theme1_hg19)
head(theme1_hg19_gene)
head(theme2_hg19)
head(theme2_hg19_gene)

head(long_short_sum)
dim(long_short_sum)

long_hg19_whole$HLAgene <- "5M"
short_hg19_whole$HLAgene <- "5M"
theme1_hg19$HLAgene <- "5M"
theme2_hg19$HLAgene <- "5M"
theme3_hg19$HLAgene <- "5M"

theme2_dedup_hg19$HLAgene <- "5M"
theme3_dedup_hg19$HLAgene <- "5M"


long_short_sum$HLAgene <- "5M"
head(long_short_sum_gene)
head(long_short_sum)
long_hg19$type <- "Longread"
short_hg19_dedup$type <- "Shortread"

long_hg19 %>% rbind(long_hg19_whole)
long_hg19 %>% left_join(long_hg19_whole) %>% count(HLAgene)
long_hg19 %>% rbind(short_hg19_dedup) %>%  #count(HLAgene)
  rbind(theme1_hg19_gene) %>% rbind(theme2_hg19_gene) %>% #head()#count(type)
  #rbind(long_short_sum_gene) %>%# head()
  select(-contig) %>%
  rbind(long_hg19_whole) %>% 
  rbind(short_hg19_whole) %>% 
  #rbind(long_short_sum) %>% #head()
  rbind(theme1_hg19) %>% rbind(theme2_hg19) %>% #head()#count(type)
  select(meandepth,numreads,ID,HLAgene,type) %>% #head()
  ggplot(aes(x=type,y=meandepth,fill=type)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 1) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") -> p1
  
long_hg19 %>% rbind(short_hg19_dedup) %>%  #count(HLAgene)
  rbind(theme1_hg19_gene) %>% rbind(theme2_hg19_gene) %>% #head()#count(type)
  #rbind(long_short_sum_gene) %>%# head()
  select(-contig) %>%
  rbind(long_hg19_whole) %>% 
  rbind(short_hg19_whole) %>% 
  #rbind(long_short_sum) %>% #head()
  rbind(theme1_hg19) %>% rbind(theme2_hg19) %>% #head()#count(type)
  select(meandepth,numreads,ID,HLAgene,type) %>% #head()
  filter(HLAgene != "5M") %>%
  ggplot(aes(x=type,y=numreads,fill=type)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 1) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") -> p2



long_hg19 %>% rbind(short_hg19_dedup) %>%  #count(HLAgene)
  rbind(theme1_hg19_gene) %>% rbind(theme2_hg19_gene) %>% #head()#count(type)
  #rbind(long_short_sum_gene) %>%# head()
  select(-contig) %>%
  rbind(long_hg19_whole) %>% 
  rbind(short_hg19_whole) %>% 
  #rbind(long_short_sum) %>% #head()
  rbind(theme1_hg19) %>% rbind(theme2_hg19) %>% #head()#count(type)
  select(meandepth,numreads,ID,HLAgene,type) %>% #head()
  filter(HLAgene != "5M") %>%
  ggplot(aes(x=type,y=numreads,fill=type)) +
  geom_boxplot() +
  facet_wrap(~HLAgene,nrow = 1) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") -> p3


p3 <- get_legend(p3)

plot_grid(
  p1,p2,p3,nrow=3,
  rel_heights = c(4, 4,0.5)
)



head(theme1_hg19_exon)
head(long_hg19_whole)

long_hg19_whole %>% rbind(short_hg19_whole) %>% 
  mutate('exon'="5M") %>% #head()
  rbind(theme1_hg19_exon %>% rbind(theme2_hg19_exon) %>% select(-contig)) %>% #head()
  filter(exon != "5M") %>%
  select(meandepth,numreads,ID,HLAgene,exon,type) %>% #head()
  ggplot(aes(x=exon,y=meandepth,fill=type)) +
  geom_boxplot() + 
  #facet_grid(exon~HLAgene)+
  #facet_grid(~HLAgene,)+
  facet_wrap(~HLAgene, ncol = 3) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")



head(theme1_hg19)
head(theme1_hg19_gene)
head(theme2_hg19)
head(theme2_hg19_gene)
head(theme2_dedup_hg19)
head(theme2_dedup_hg19_gene)
head(theme3_hg19)
head(theme3_hg19_gene)
head(theme3_dedup_hg19)
head(theme3_dedup_hg19_gene)

theme1_hg19_gene %>% rbind(theme2_hg19_gene) %>% rbind(theme2_dedup_hg19_gene) %>%
  rbind(theme3_hg19_gene) %>% rbind(theme3_dedup_hg19_gene) %>% 
  select(-contig) %>%
  rbind(theme1_hg19) %>% rbind(theme2_hg19) %>% rbind(theme3_hg19) %>%
  rbind(theme2_dedup_hg19) %>% rbind(theme3_dedup_hg19) %>% #head()#count(type,HLAgene) %>% head()
  select(ID,type,HLAgene,meandepth) -> merge_hg19_all

head(merge_hg19_all)  

merge_hg19_all %>% 
  #select(meandepth,numreads,ID,HLAgene,exon,type) %>% #head()
  ggplot(aes(x=type,y=meandepth,fill=type)) +
  geom_boxplot() + 
  #facet_grid(exon~HLAgene)+
  #facet_grid(~HLAgene,)+
  facet_wrap(~HLAgene, ncol = 9) + 
    theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")


