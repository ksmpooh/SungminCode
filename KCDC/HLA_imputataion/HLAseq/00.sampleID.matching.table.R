### 일반용역 데이터 아이디 정리
library(readxl)
library(stringr)
library(tidyverse)
library(stringr)

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
#table(hlaseq$cell)
table(hlaseq$`Bio Sample Name`)

#MHC_Illumina_Rawdata
#MHC_PacBio_Rawdata

## long read
#./2nd_Cell_CCS/2nd_Cell.bc1001--bc1001.bam
hlaseq$cell
hlaseq$Cell[1:12]
#hlaseq[(2:12),]$Cell <- 1
hlaseq$cell <- hlaseq$Cell[1]
hlaseq$cell[13:24] <- hlaseq$Cell[13]
hlaseq$cell[25:36] <- hlaseq$Cell[25]
hlaseq$cell[37:48] <- hlaseq$Cell[37]
hlaseq$cell[49:60] <- hlaseq$Cell[49]

head(hlaseq)


str_split_fixed(hlaseq$cell," ",3)[,2] %>% head()
str_split_fixed(hlaseq$cell," ",3) %>% head()


table(hlaseq$cell)
head(hlaseq)
df_sub <- df %>% select(KBA_ID.2020,HLAseqID) %>% mutate(index = ifelse(grepl("2020HLA",HLAseqID),HLAseqID,str_sub(KBA_ID.2020,-4)))
hlaseq %>% select(`Bio Sample Name`,cell,`Barcode Name`) %>% rename(HLA_pro_ID = `Bio Sample Name`,barcode_Name = `Barcode Name`) %>% 
  merge(df_sub,by.x = "HLA_pro_ID",by.y="index") %>% #head()
  mutate(Longread_filename = paste0(str_split_fixed(cell," ",3)[,1],"_",str_split_fixed(cell," ",3)[,2],".",barcode_Name)) %>% 
  mutate(Longread_filePath = paste0("./MHC_PacBio_Rawdata/",str_split_fixed(cell," ",3)[,1],"_",str_split_fixed(cell," ",3)[,2],"_CCS/",Longread_filename)) %>%
  show()



## short read
#1003_S71_L002_R1_001.fastq.gz
