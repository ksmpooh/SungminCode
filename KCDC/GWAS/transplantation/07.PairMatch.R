setwd("c:/Users/user/Desktop/KCDC/FinalKing_JG/")
df1 <-read.table("04.newKing.6772/01.allsnp/02.related.2ndDegree/related.2ndDegree.Result.kin0",header = T)
df2 <-read.table("04.newKing.6772/02.filter/02.related.2ndDegree/filter.2ndDegree.kin0",header = T)
df3 <- read.table("04.newKing.6772/03.pruning/02.related.2ndDegree/JG.pruning50.5.0.1.2ndDegree.kin0",header = T)

table(df1$InfType)
df1 <- df1[df1$InfType %in% c('2nd','Dup/MZ','FS','PO'),]
df2 <- df2[df2$InfType %in% c('2nd','Dup/MZ','FS','PO'),]
df3 <- df3[df3$InfType %in% c('2nd','Dup/MZ','FS','PO'),]
library(readr)
library(dplyr)
real <- readr::read_csv("../transplantation/final_sample_info/pairtable/KOTRY_KCHIP_ID_full_20200306.csv",col_names = T)
real <- as.data.frame(real)

head(real)

real_pair <- real[,c(1,8,9,10,11,12,19,20,21,22)]

head(real_pair)
colnames(real_pair)

real_pair <- real_pair[,c(1,4,6,10)]
real_pair <- real_pair[,c(1,6)]

real_pair <-na.omit(real_pair)
head(real_pair)



head(df1)
df1 <-df1[,c(1,3,14)]
df2 <-df2[,c(1,3,14)]
df3 <-df3[,c(1,3,14)]

colnames(df1)[1:2]<-colnames(real_pair)
colnames(df2)[1:2]<-colnames(real_pair)
colnames(df3)[1:2]<-colnames(real_pair)
head(df1)
head(df2)
head(df3)


ref <- merge(real_pair,df1,by=c("KCHIP_ID","KCHIP_ID1"))
ref <- merge(real_pair,df2,by=c("KCHIP_ID","KCHIP_ID1"))
ref <- merge(real_pair,df3,by=c("KCHIP_ID","KCHIP_ID1"))
table(ref$InfType)
head(ref)
tail(ref)




###  연구노트 용

# 예측결과 읽기
king <- read.table("king_result.kin0",header = T)
# pair와 예측관계(Inftype) subset
king <- king[,c("ID1","ID2","InfType")] 
# real pair 읽기
real <- read.talbe("real_pair.txt",header = T)

# real의 ID와 king ID를 동일하게 설정
colnames(king)[1:2] <- colnames(real)

# merge 하고 pair 숫자 확인
df <- merge(real,king,by=c(colnames(king)[1:2)))
nrow(df)
table(df$InfType)






