setwd("C:/Users/user/Desktop/KCDC/transplantation/")
library(readr)
library(dplyr)
real <- readr::read_csv("final_sample_info/pairtable/KOTRY_KCHIP_ID_full_20200306.csv",col_names = T)
real <- as.data.frame(real)

head(real)

real_pair <- real[,c(1,8,9,10,11,12,19,20,21,22)]

head(real_pair)
colnames(real_pair)

real_pair <- real_pair[,c(1,6,4)]
#real_pair <- real_pair[,c(1,6)]

real_pair <-na.omit(real_pair)
head(real_pair)
table(real_pair$Rela_Pair)
real_pair$Rela_Pair<-factor(real_pair$Rela_Pair,labels = c("PO","PO","Dup/MZ","FS","Spouse","Relative","unrelated","Dual graft"))
write.csv(real_pair,"../FinalKing_JG/04.newKing.6772/00.ref/KOTRY_KCHIP_Realpair.onlyPairInfo.20201209.csv",row.names = F,quote = F)
