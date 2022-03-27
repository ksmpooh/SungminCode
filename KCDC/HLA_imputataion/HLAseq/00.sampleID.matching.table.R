### 일반용역 데이터 아이디 정리
library(readxl)
library(stringr)
library(tidyverse)

setwd("~/Desktop/KCDC/일반용역/")



df <- read.table("HLA.typing.Final.result.20211126.xlsx")
df <- read_excel("HLAtyping.265pairtable_with(QC.typing)_20211028.xls",sheet = 1)
head(df)

hlaseq <- read_excel("2021_HLAseq_neurogen/HLA_seq_2021/long/20220125.MHC.HLA.1-5.Summary.xlsx",sheet = 1)
hlaseq_short <- read.table("2021_HLAseq_neurogen/HLA_seq_2021/short/file.list.txt",header=F)
hlaseq <- hlaseq[hlaseq$`Bio Sample Name` != "No Name",]
head(hlaseq)
head(hlaseq_short)
table(hlaseq$Cell)
table(hlaseq$cell)
table(hlaseq$`Bio Sample Name`)


## long read
hlaseq$Cell[1:12]
#hlaseq[(2:12),]$Cell <- 1
hlaseq$cell <- hlaseq$Cell[1]
hlaseq$cell[13:24] <- hlaseq$Cell[13]
hlaseq$cell[25:36] <- hlaseq$Cell[25]
hlaseq$cell[37:48] <- hlaseq$Cell[37]
hlaseq$cell[49:60] <- hlaseq$Cell[49]


head(hlaseq)

table(hlaseq$cell)
hlaseq %>% select(`Bio Sample Name`,cell,`Barcode Name`) %>% head()


## short read
