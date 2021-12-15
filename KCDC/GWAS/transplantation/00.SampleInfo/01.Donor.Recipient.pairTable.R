### Data reference
setwd("~/Desktop/KCDC/일반용역/")
library(readxl)
library(writexl)
library(tidyverse)
library(stringr)
ref <- read_excel("HLAtyping.265pairtable_with(QC.typing)_20211028.xls",sheet = 1) %>%
  select(KBA_ID.2019,oriID.2019,HLAID.2019,KBA_ID.2020,oriID.2020,HLAID.2020) %>% as.data.frame() 
colnames(ref)                  
head(ref)
ref$HLAID.2019 <- str_replace_all(ref$HLAID.2019,"H","CDC")

#####HLA typing


setwd("~/Desktop/KCDC/transplantation/00.sampleInfo/")
dis <- read.table("2019/last.sample.info_20201207.txt",header = T) %>%
  select(KBA_ID,OriID,type) %>% mutate(Prod = "2019") %>% mutate(ref = str_sub(OriID,3,-1)) %>%
  filter(type == "KD" | type == "KR")
head(dis)
table(dis$type)

rep <- read.csv("2020/2020_JG_KCHIP_ID.csv") %>%
  select(KID,ID,Type) %>% rename(KBA_ID = KID,OriID = ID,type = Type) %>%
  mutate(Prod = "2020") %>% mutate(ref = str_sub(OriID,3,-1)) %>%
  filter(type == "KD" | type == "KR")
head(rep)
table(rep$type)

table(dis$type) + table(rep$type)


KD <- rbind(dis %>% filter(type == "KD"),rep %>% filter(type=="KD"))
KR <- rbind(dis %>% filter(type == "KR"),rep %>% filter(type=="KR"))
head(KD)


df <- merge(KR,KD,by = 'ref',all = T) %>% mutate(Check = ifelse(Prod.x == Prod.y,1,0))
table(df$Check)
df
