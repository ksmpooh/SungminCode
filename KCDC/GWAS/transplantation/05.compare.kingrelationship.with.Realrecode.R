library(readr)
library(dplyr)
setwd("c:/Users/user/Desktop/KCDC/transplantation/")
real <- readr::read_csv("final_sample_info/pairtable/KOTRY_KCHIP_ID_full_20200306.csv",col_names = T)
real <- as.data.frame(real)
head(real)

table(real$Rela_Pair)



real_pair <- real[,c(1,8,9,10,11,12,19,20,21,22)]

head(real_pair)
colnames(real_pair)

#real_pair <- real_pair[,c(1,4,6,10)]
real_pair <- real_pair[,c(1,6)]

real_pair <-na.omit(real_pair)
#real_pair <-real_pair[,]
result <- read.table("2nd/all.king/2ndkingResult.withType.txt",header = T)
head(result)
############2020.03.20
table(result$FID1.type)
result <- result[!result$FID1.type == 'control',]
result[result$FID1 == "NIH19KT0017",]
result[result$FID2 == "NIH19KT0017",]
#####
result_pair <- result[,c(1,3)]
head(result_pair)

head(real_pair)
colnames(real_pair) <- c("FID1","FID2")
#colnames(real_pair)[1] <- "FID1"
#colnames(real_pair)[3] <- "FID2"
################20200320
result_pair$FID1 <- as.character(result_pair$FID1)
result_pair$FID2 <- as.character(result_pair$FID2)
#real_pair$FID1 <- as.factor(real_pair$FID1)
#real_pair$FID2 <- as.factor(real_pair$FID2)
###################
head(result_pair)

real_pair[real_pair$FID1 == "NIH19KT0017",]
str(real_pair)
str(result_pair)


df[df$FID1 == "NIH19KT0017",]
real_pair[real_pair$FID1 == "NIH19KT0017",]
result_pair[result_pair$FID1 == "NIH19KT0017",]

head(result_pair)
head(real_pair)
df <- real_pair %>%anti_join(result_pair,by=c("FID1","FID2"))
head(df)
colnames(df) <- c("FID2","FID1")
#colnames(real_pair)[1] <- "FID2"
#colnames(real_pair)[3] <- "FID1"

df <- df %>%anti_join(result_pair,by=c("FID2","FID1"))
head(df)
colnames(df) <- c("FID1","FID2")
df[df$FID1 == "NIH19KT0017",]
head(df)

colnames(df) <- c("KCHIP_ID","KCHIP_ID1")
head(real)



rmlist <- read.table("final_sample_info/last.sample.info.txt",header = T)

head(rmlist)
table(rmlist$state)
rmlist <- rmlist[rmlist$state == "remove",]
colnames(rmlist)[1]<-"KCHIP_ID"

df <- df %>% anti_join(rmlist,by = "KCHIP_ID")

colnames(rmlist)[1]<-"KCHIP_ID1"

df <- df %>% anti_join(rmlist,by = "KCHIP_ID1")
df[df$KCHIP_ID == "NIH19KT0017",]


df <- merge(df,real,by = "KCHIP_ID")
df <- df[,c(-2)]
colnames(df)[12] <- "KCHIP_ID1"
head(df)
df[df$KCHIP_ID == "NIH19KT0017",]

write.table(df,"final_sample_info/pairtable/not.predicted.by.king.txt",col.names = T,row.names = F,quote = F)
head(df)
table(df$Rela_Pair)
table(df$Rela_Pair1)

table(real$Rela_Pair)
df[df$Rela_Pair == 1,]$KCHIP_ID




a<-df$KCHIP_ID
b <- factor(df$KCHIP_ID1)
c <- rbind(a,b)
head(c)


wr