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
hlaseq_short <- read.table("2021_HLAseq_neurogen/HLA_seq_2021/short/R1.file.list.txt",header=F)
hlaseq_short2 <- read.table("2021_HLAseq_neurogen/HLA_seq_2021/short/R2.file.list.txt",header=F)
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

str_sub(df$KBA_ID.2020,-4,-3)
table(hlaseq$cell)
head(hlaseq)
df_sub <- df %>% select(KBA_ID.2020,HLAseqID) %>% 
  mutate(index = ifelse(grepl("2020HLA",HLAseqID),HLAseqID,ifelse(str_sub(KBA_ID.2020,-4,-4) == "0",str_sub(KBA_ID.2020,-3),str_sub(KBA_ID.2020,-4))))
longread <- hlaseq %>% select(`Bio Sample Name`,cell,`Barcode Name`) %>% rename(longread_ID = `Bio Sample Name`,barcode_Name = `Barcode Name`) %>% 
  merge(df_sub,by.x = "longread_ID",by.y="index",all.x=T) %>% 
  mutate(Longread_filename = paste0(str_split_fixed(cell," ",3)[,1],"_",str_split_fixed(cell," ",3)[,2],".",barcode_Name)) %>% 
  mutate(Longread_filePath = paste0("./MHC_PacBio_Rawdata/",str_split_fixed(cell," ",3)[,1],"_",str_split_fixed(cell," ",3)[,2],"_CCS/",Longread_filename)) #%>% 
  #head()
  
  


head(df)
df_index <- df %>% select(KBA_ID.2020) %>% 
  mutate(index = ifelse(str_sub(KBA_ID.2020,-4,-4) == "0",str_sub(KBA_ID.2020,-3),str_sub(KBA_ID.2020,-4)))
## short read
#1003_S71_L002_R1_001.fastq.gz
head(hlaseq_short)
hlaseq_short2 <- hlaseq_short2 %>%rename(Shortread_filename_R2 = V1) %>% mutate(index = str_split_fixed(Shortread_filename_R2,"_",5)[,1])
shortread <- hlaseq_short %>%rename(Shortread_filename_R1 = V1) %>% mutate(index = str_split_fixed(Shortread_filename_R1,"_",5)[,1]) %>%
  merge(hlaseq_short2,by='index') %>%
  merge(df_index,by='index') %>% rename(KBA_ID = KBA_ID.2020,shortread_ID = index) %>% 
  select(shortread_ID,KBA_ID,Shortread_filename_R1,Shortread_filename_R2)

#write.table(shortread,"2021_HLAseq_neurogen/HLA_seq_2021/short/Shortread.filename.withNIHID.txt",col.names = T,row.names = F,quote = F,sep = "\t")


head(longread)
head(shortread)

merge(shortread,longread %>% select(longread_ID,KBA_ID.2020,Longread_filename,Longread_filePath),by.x="KBA_ID",by.y = "KBA_ID.2020") %>% 
  head()
