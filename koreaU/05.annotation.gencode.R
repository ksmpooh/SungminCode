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
BiocManager::repositories('CRAN: https://cran.rstudio.com/')
BiocManager::install("rtracklayer")
require(BiocManager)
install.packages('installr')
require(installr)
updateR()



d = import('~/Desktop/KU/2021_Fall/biovisual/gencode.v32.basic.annotation.gtf')


d = as.data.frame(d)

d %>% filter(type='transcript') %>% nrow()

d %>% count(strand)
d %>% filter(seqnames == 'chr3', type == 'transcript') %>% nrow()