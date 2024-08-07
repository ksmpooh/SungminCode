### imputation 결과에 HLA typing (NGS)결과 변경
library(tidyverse)
library(stringr)

setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/")

ref <-read.table("HLA.typing_2digit.txt",header = T,sep = "\t",colClasses = c('character'))
ref <-read.table("HLA.typing_4digit.txt",header = T,sep = "\t",colClasses = c('character'))

head(ref)

setwd("~/Desktop/KCDC/HLAimputation/03.HLA_matching/")

#df <- read.csv("2019_Kidney_HLAmathing_table_2digit.csv")
df <- read.table("2019_Kidney_HLAmathing_table_2digit.txt",header = T,sep = "\t",colClasses = c('character'))
df <- read.table("2019_Kidney_HLAmathing_table_4digit.txt",header = T,sep = "\t",colClasses = c('character'))
head(df)



####
colnames(df)
R <- ref %>% filter(ref$KBA_ID %in% df$KR_KBA_ID) %>% mutate(ref = substr(OriID,3,7)) %>% select(ref,KBA_ID,OriID,A.1,A.2,B.1,B.2,C.1,C.2,DRB1.1,DRB1.2,DQA1.1,DQA1.2,DQB1.1,DQB1.2,DPA1.1,DPA1.2,DPB1.1,DPB1.2)
D <- ref %>% filter(ref$KBA_ID %in% df$KD_KBA_ID) %>% mutate(ref = substr(OriID,3,7)) %>% select(ref,KBA_ID,OriID,A.1,A.2,B.1,B.2,C.1,C.2,DRB1.1,DRB1.2,DQA1.1,DQA1.2,DQB1.1,DQB1.2,DPA1.1,DPA1.2,DPB1.1,DPB1.2)

organ = "K"
#paste0(organ,"R_A.1")
colnames(R)[4:19] <- c(paste0(organ,"R_A.1"),paste0(organ,"R_A.2"),paste0(organ,"R_B.1"),paste0(organ,"R_B.2"),paste0(organ,"R_C.1"),paste0(organ,"R_C.2"),paste0(organ,"R_DRB1.1"),paste0(organ,"R_DRB1.2"),paste0(organ,"R_DPA1.1"),paste0(organ,"R_DPA1.2"),paste0(organ,"R_DPB1.1"),paste0(organ,"R_DPB1.2"),paste0(organ,"R_DQA1.1"),paste0(organ,"R_DQA1.2"),paste0(organ,"R_DQB1.1"),paste0(organ,"R_DQB1.2"))
colnames(D)[4:19] <- c(paste0(organ,"D_A.1"),paste0(organ,"D_A.2"),paste0(organ,"D_B.1"),paste0(organ,"D_B.2"),paste0(organ,"D_C.1"),paste0(organ,"D_C.2"),paste0(organ,"D_DRB1.1"),paste0(organ,"D_DRB1.2"),paste0(organ,"D_DPA1.1"),paste0(organ,"D_DPA1.2"),paste0(organ,"D_DPB1.1"),paste0(organ,"D_DPB1.2"),paste0(organ,"D_DQA1.1"),paste0(organ,"D_DQA1.2"),paste0(organ,"D_DQB1.1"),paste0(organ,"D_DQB1.2"))

head(R)
head(D)
head(df)
out <- merge(R %>% rename(KR_KBA_ID = 'KBA_ID',KR_OriID = 'OriID'),
             D %>% rename(KD_KBA_ID = 'KBA_ID',KD_OriID = 'OriID'),by='ref') 

colnames(out)
out <- out[,c(1:3,20:21,4:19,22:37)]
head(out)

#ori <- read.table("~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_2digit.txt",header = T,sep = "\t")
#ori <- read.table("~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_4digit.txt",header = T,sep = "\t")

head(df)
df <- df %>% filter(!(df$KR_KBA_ID %in% out$KR_KBA_ID)) %>% rbind(out)
df %>% filter(!(df$KR_KBA_ID %in% out$KR_KBA_ID)) %>% rbind(out)

write.table(df,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_2digit_changeNGStyping.txt",row.names = F,quote = F,sep = "\t")
write.table(df,"~/Desktop/KCDC/HLAimputation/03.HLA_matching/2019_Kidney_HLAmathing_table_4digit_changeNGStyping.txt",row.names = F,quote = F,sep = "\t")






head(out)
write.table(out %>% select(KR_KBA_ID,KD_KBA_ID),"~/Desktop/KCDC/HLAimputation/03.HLA_matching/NGS.HLAtyping.ID.list.txt",row.names = F,quote = F,sep = "\t")


