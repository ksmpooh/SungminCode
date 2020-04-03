setwd("C:/Users/user/Desktop/KCDC/transplantation/")
library(readr)
library(dplyr)
real <- readr::read_csv("final_sample_info/pairtable/KOTRY_KCHIP_ID_full_20200306.csv",col_names = T)
real <- as.data.frame(real)

head(real)

table(real$Rela_Pair)
real_pair <- real[,c(1,8,9,10,11,12,19,20,21,22)]

head(real_pair)
colnames(real_pair)

real_pair <- real_pair[,c(1,4,6,10)]
#real_pair <- real_pair[,c(1,6)]

real_pair <-na.omit(real_pair)


###04에서 작업해야함..작업 후 나오는 fianl result
result <- read.csv("king/Final/Final.kingResult.allsample.2ndDegree.allsnp.withType.csv",header = T)
result <- result[!result$FID1.type == 'control',]
head(result)
length(result$InfType)
#####
result_pair <- result[,c(1,2,4,5,7)]
#result_pair <- result[,c(1,4)]
head(result_pair)
head(real_pair)
#colnames(real_pair) <- c("FID1","FID2")
result_pair$FID1 <- as.character(result_pair$FID1)
result_pair$FID2 <- as.character(result_pair$FID2)


df <- merge(real_pair,result_pair,by.x = "KCHIP_ID",by.y = "FID1",all.y = T)
head(df)
df1 <- merge(real_pair,result_pair,by.x = "KCHIP_ID1",by.y = "FID2",all.y = T)
df2 <- merge(df1,df)
head(df2)
#head(df1)
#df <- real_pair %>%inner_join(result_pair,by=c("FID1","FID2"))
#colnames(real_pair)<-c("FID2","FID1")
#df1 <- real_pair %>%inner_join(result_pair,by=c("FID1","FID2"))
#df <- real_pair %>% anti_join(result_pair,by=c("FID1","FID2"))
#df <- real_pair %>%inner_join(result_pair,by=c("FID1","FID2"))

#colnames(df) <- c("FID2","FID1")
#colnames(real_pair)[1] <- "FID2"
#colnames(real_pair)[3] <- "FID1"

write.csv(df,)
