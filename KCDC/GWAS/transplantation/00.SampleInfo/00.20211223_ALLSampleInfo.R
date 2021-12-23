library(readxl)
library(stringr)
library(tidyverse)
library(writexl)
df <- read_excel("~/Desktop/KCDC/transplantation/00.sampleInfo/QCed.list/KOTRY_samples_20211114_checked.xlsx",sheet = 1)
head(df)

df$type.1 <- substr(df$...3,1,2)
df$type.2 <- substr(df$...7,1,2)

df%>% filter(df$수여자 == "2019") %>% count(type.1)
df%>% filter(df$수여자 == "2019") %>% count(type.2)


table(df$수여자)
table(df$공여자)

table(df$type.1)
table(df$type.2)


setwd("~/Desktop/KCDC/transplantation/00.sampleInfo/")
ref <-read_excel("~/Desktop/KCDC/transplantation/00.sampleInfo/ALL(2019and2020).sampleID.xlsx",sheet = 1)
sum_ref <- read.table("2020/JG.2020.NIH.ID.txt")
head(ref)
head(sum_ref)

dis <- read.table("2019/last.sample.info_20201207.txt",header = T) %>% select(KBA_ID,TubeID) %>% rename(bCODE = "TubeID")
head(dis)
rep <- read.csv("2020/2020_JG_KCHIP_ID.csv",header = T) %>% select(KID,bCODE) %>% rename(KBA_ID = "KID") %>% filter(KBA_ID %in% sum_ref$V1)
head(rep)


bcode <- rbind(dis,rep)
head(bcode)

ref <- merge(ref,bcode)
ref$QC <- "QCout"

KR_2019 <- read.table("QCed.list/KR.2019.QCed.list.txt")
KR_2020 <- read.table("QCed.list/KR.rep.2020.QCed.list.txt")
KD_2019 <- read.table("QCed.list/KD.2019.QCed.list.txt")
KD_2020 <- read.table("QCed.list/2020.KD.QCed.ID.txt")
LR_2019 <- read.table("QCed.list/LR.2019.QCed.list.txt")
head(LR_2019)
head(KR_2019)
head(KR_2020)
head(KD_2019)
head(KD_2020)

QCed <- rbind(KR_2019,KR_2020,KD_2019,KD_2020,LR_2019)

ref[ref$KBA_ID %in% QCed$V1,]$QC <- "QCin"
ref[ref$type == "LD",]$QC <- "nonQC"
ref[ref$type == "LR" & ref$Prod == "2020",]$QC <- "nonQC"
head(ref)
table(ref$QC)

#write_xlsx(ref,"QClist/KOTRY_KCHIPprod_ALLsample.SampleInfo.20211223.xlsx")

rec <- ref %>% filter(type == 'LR' | type == 'KR')
don <- ref %>% filter(type == 'LD' | type == 'KD')

out <- merge(rec,don,by='ref',all = T)
head(out)



write_xlsx(out,"QClist/KOTRY_KCHIPprod_ALLsample.SampleInfo_pairtable.20211223.xlsx")

table(ref$QC)
sum(table(ref$QC))


