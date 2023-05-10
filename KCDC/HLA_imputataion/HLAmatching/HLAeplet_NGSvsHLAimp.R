### HLA eplet about NGS vs
### 1. HLA NGS typing check
### 2. HLA typing vs HLA imputation


##### 20230508 HLA imp result vs NGS result
ngs <- read.table("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_529sample_IMGT3320convert.txt",header = T)
ngs %>% select(-Sample) -> ngs
head(ngs)
head(dr)


dr %>% merge(ngs,by.x = 'KR',by.y = "ID") %>%
  merge(ngs,by.x = 'KD',by.y = "ID") -> ngs_out

head(ngs_out)
tcol <- grep("NGS_",colnames(ngs_out))
tcol
ngs_out[,tcol] <- lapply(ngs_out[tcol], function(x)
  paste0(str_split_fixed(x,":",3)[,1],":",str_split_fixed(x,":",3)[,2])
)
head(ngs_out)
ngs_out$x <- "x"
ngs_out$x1 <- "x"
ngs_out$x2 <- "x"
ngs_out$x3 <- "x"

ngs_out %>%
  select(KR,NGS_A.1.x,NGS_A.2.x,NGS_B.1.x,NGS_B.2.x,NGS_C.1.x,NGS_C.2.x,KD,NGS_A.1.y,NGS_A.2.y,NGS_B.1.y,NGS_B.2.y,NGS_C.1.y,NGS_C.2.y) -> c1
#writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class1.xlsx",col_names = F,) -> c1
ngs_out %>%
  select(KR,NGS_DRB1.1.x,NGS_DRB1.2.x,x,x1,NGS_DQB1.1.x,NGS_DQB1.2.x,NGS_DQA1.1.x,NGS_DQA1.2.x,NGS_DPB1.1.x,NGS_DPB1.2.x,NGS_DPA1.1.x,NGS_DPA1.2.x,KD,NGS_DRB1.1.y,NGS_DRB1.2.y,x2,x3,NGS_DQB1.1.y,NGS_DQB1.2.y,NGS_DQA1.1.y,NGS_DQA1.2.y,NGS_DPB1.1.y,NGS_DPB1.2.y,NGS_DPA1.1.y,NGS_DPA1.2.y) -> c2
#writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class2.xlsx",col_names = F) -> c2
dataset_names <- list("classI" = c1, "classII" = c2)
writexl::write_xlsx(dataset_names,"test_mtsheet.xlsx")




library(readxlsb)
setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

df <- read_xlsb("NGS_eplet/ABC_Eplet_Matching_4.0_HLAtyping.xlsb",sheet = "Check HLA")
head(df)

colnames(df)
select(grep("Check",colnames(df)))
df %>% select(grep("Check",colnames(df))) %>% #head()
  filter(!row_number() %in% c(1,2)) %>% #head()
  pivot_longer(1:12,names_to = "check",values_to = "type") %>% #head()
  na.omit() %>% select(type) %>% unique() -> c1

df <- read_xlsb("NGS_eplet/DRDQDP_Eplet_Matching_3.1_HLAtyping.xlsb",sheet = "Check HLA")
#colnames(df)
df %>% select(grep("check",colnames(df))) %>% #head()
  filter(!row_number() %in% c(1,2)) %>% #head()
  pivot_longer(1:24,names_to = "check",values_to = "type") %>% #head()
  na.omit() %>% select(type) %>% #head()
  unique() -> c2


head(c1)
head(c2)

c1$class <- "classI"
c2$class <- "classII"

c <- rbind(c1,c2)
head(c)

c %>% select(class,type) %>% arrange(class,type) %>% na.omit() %>% #head()
  write_xlsx("NGS_eplet/Check.HLAtypelist.foreplet.xlsx")




### 2. HLA typing vs HLA imputation

library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")
