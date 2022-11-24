library(readxl)
library(tidyverse)
library(stringr)
setwd("~/Desktop/KCDC/KCHIP_open/임신합병증/")

df <- read.table("check.MZ.sample.imiss",header = T)
head(df)
#str_replace_all(df$FID,".CEL","")
df$ID <- str_split_fixed(str_replace_all(df$FID,".CEL",""),"_",6)[,6]

ref1 <- read_excel("@preg_MZforQC.check.xlsx",sheet = 1)
ref2 <- read_excel("@preg_MZforQC.check.xlsx",sheet = 2)
head(ref1)
head(ref2)
rm1 <- ref1 %>% rbind(ref2) %>% select(ID1,`바이오뱅크과 및 역학과 확인`,후속조치) %>% rename("ID" = ID1)
rm2 <- ref1 %>% rbind(ref2) %>% select(ID2,`바이오뱅크과 및 역학과 확인`,후속조치) %>% rename("ID" = ID2)
rmall <- rbind(rm1,rm2)
rmall
head(df)

ref3 <- read_excel("@preg_MZforQC.check.xlsx",sheet = 3)
ref4 <- read_excel("@preg_MZforQC.check.xlsx",sheet = 4)
ref3 %>% rbind(ref4) %>% select(ID1,ID2,`바이오뱅크과 및 역학과 확인`,후속조치) %>% #rename("ID" = ID1) %>%
  left_join(df %>% select(ID,F_MISS) %>% rename("ID1" = ID),by='ID1') %>% rename("ID1_MISS" = F_MISS) %>%
  left_join(df %>% select(ID,F_MISS) %>% rename("ID2" = ID),by='ID2') %>% rename("ID2_MISS" = F_MISS) %>% #mutate(ID = ifelse(ID1_MISS >= ID2_MISS,ID2,ID1)) ->t
  mutate(ID = ifelse(ID1_MISS >= ID2_MISS,ID2,ID1)) %>% select(ID,`바이오뱅크과 및 역학과 확인`,후속조치) -> rmselect

head(rmselect)
head(df)
rmall %>% rbind(rmselect) %>% na.omit() %>% inner_join(df %>% select(ID,FID),by='ID') -> rmlist
rmlist
write.table(na.omit(rmlist),"rmlist.MZlistafterBioBankcheck_20221124.txt",col.names = T,row.names = F,quote = F,sep = "\t")
writexl::write_xlsx(na.omit(rmlist),"rmlist.MZlistafterBioBankcheck_20221124.xlsx")

na.omit(rmlist) %>% count(`바이오뱅크과 및 역학과 확인`)
