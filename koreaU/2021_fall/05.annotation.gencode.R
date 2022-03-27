## 20211011 수업

setwd("~/Desktop/KU/2021_Fall/biovisual/")

library(tidyverse)
d = read_delim('gencode.v32.basic.annotation.gtf',
               delim='\t', skip = 5, progress = F, 
               col_names = F)
# read.csv, read.delim 보다 빠르다...ㄷ

head(d)
str(d)


install.packages("BiocManager")
library(BiocManager)
BiocManager::install("rtracklayer")

library(rtracklayer)
library(tidyverse)
d = import('~/Desktop/KU/2021_Fall/biovisual/gencode.v32.basic.annotation.gtf')


d = as.data.frame(d)
cols = c('chrom', 'source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 'info')
head(d)
d %>% filter(type=="transcript") %>% nrow()

d %>% count(strand)
d %>% filter(seqnames == 'chr3', type == 'transcript') %>% nrow()

d %>% filter(type == 'exon') %>% group_by(gene_type) %>% count()

head(d)

## Define TTS upstream promoter

d %>% filter(type == 'transcript') %>% mutate(position = start - 2000) %>% head()
d %>% filter(type == 'transcript') %>%
  mutate(TTS_position = ifelse(strand == "+", start,end), 
         promoter_start1 = ifelse(strand == "+",TTS_position - 2000,TTS_position + 2000),
         promoter_end1 = TTS_position,
         promoter_start = ifelse(promoter_start1 < promoter_end1,promoter_start1,promoter_end1),
         promoter_end = ifelse(strand == '+',TTS_position,promoter_start1)) %>%
  select(gene_name,transcript_id, seqnames,TTS_position,strand,promoter_start,promoter_end) %>% head()


## 20211018

library(tidyverse)
d = read_delim('~/Desktop/KU/2021_Fall/biovisual/gencode.v32.basic.annotation.gtf',
               delim='\t', skip = 5, progress = F, 
               col_names = F)

d = as.data.frame(d)
head(d)
cols = c('chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'info')

d = read_delim('~/Desktop/KU/2021_Fall/biovisual/gencode.v32.basic.annotation.gtf', 
               delim='\t', skip = 5, 
               progress = F,
               col_names = cols)
d = as.data.frame(d)

d %>% filter(type == 'UTR')
d %>% filter(type == 'UTR')%>%
  group_by(chrom = seqnomaes) %>%
  summarize(mean_legnth = ment(width)) %>%
  ggplot(aes(chrom, mean_length)) + geom_bar(stat = 'identity') + coord_flip()


d %>% filter(type== 'transcript') %>%
  group_by(chrom=seqnames) %>%
  summarize(mean_length = mean(width))

## 염색체 별로 유전자의 갯수를 찾으라
head(d)
d %>% filter(type == 'gene') %>%
  group_by(chrom=seqnames) %>% count(gene_name) %>%
  arrange(desc())

d %>% filter(type == 'gene') %>% head()

## 염색체 별로 유전자의 갯수를 찾으라 액손을 많이 갖는 유전자 상휘 10개를 살펴보자

d %>% filter(type == 'exon') %>%
  group_by(seqnames) %>%
  count(gene_name) %>%
  arrange(des(n))


d %>% select(chrom,type)



table(d$type)

d %>% 


