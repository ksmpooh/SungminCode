##### 20200811 JG.1st.SNPOlist.king result compare with clinical data


setwd("c:/Users/user/Desktop/KCDC/FinalKing_JG/04.newKing.6772/")
#ref<-read.csv("00.ref/KOTRY_KCHIP_ID_full_20200306.csv",header = T)


#### pair table 
library(readr)
library(dplyr)
real <- readr::read_csv("00.ref/KOTRY_KCHIP_ID_full_20200306.csv",col_names = T)
real <- as.data.frame(real)

head(real)

table(real$Rela_Pair)
real_pair <- real[,c(1,8,9,10,11,12,19,20,21,22)]

head(real_pair)
colnames(real_pair)

real_pair <- real_pair[,c(1,6,10)]
#real_pair <- real_pair[,c(1,6)]
real_pair <-na.omit(real_pair)
table(real$Rela_Pair)

#real_pair[real_pair$Rela_Pair == '3',]
#### pair table new
#real <- read.csv("../01.1stQCed/KOTRY_KCHIP_ID_full_20200403.csv")

##### king result table with type (04_1.R)

#king <- read.table("01.allsnp/02.related.2ndDegree/related.2ndDegree.Result.kin0",header = T)
#king <- read.table("02.filter/02.related.2ndDegree/filter.2ndDegree.kin0",header = T)
king <- read.table("03.pruning/02.related.2ndDegree/JG.pruning50.5.0.1.2ndDegree.kin0",header = T)
sample.info <- read.table("../../transplantation/final_sample_info/last.sample.info.txt",header = T)
head(sample.info)
sample.info <- subset(sample.info,select = c("NewID","OriID"))


LR <- sample.info[grep("^LR",sample.info$OriID),]
KR <- sample.info[grep("^KR",sample.info$OriID),]
LD <- sample.info[grep("^LD",sample.info$OriID),]
KD <- sample.info[grep("^KD",sample.info$OriID),]

LR$type <- 'LR'
KR$type <- 'KR'
LD$type <- 'LD'
KD$type <- 'KD'

sample.info <- rbind(LR,KR)
sample.info <- rbind(sample.info,KD)
sample.info <- rbind(sample.info,LD)
head(sample.info)
sample.info <- subset(sample.info,select = c("NewID","type"))

head(king)
king <- subset(king,select = c("FID1","FID2","InfType"))

df <- merge(king,sample.info,by.y = "NewID",by.x = "FID1",all.x = T)
head(df)
colnames(df)[4] <- "FID1.type"
df <- merge(df,sample.info,by.y = "NewID",by.x = "FID2",all.x = T)
colnames(df)[5] <- "FID2.type"

df$FID1.type <- as.factor(df$FID1.type)
df$FID2.type <- as.factor(df$FID2.type)
str(df)
head(df)
tail(df)
levels(df$FID1.type)

#levels <- levels(df$FID1.type)

#levels[length(levels) + 1] <- "control"

# refactor Species to include "None" as a factor level
# and replace NA with "None"
#df$FID1.type <- factor(df$FID1.type, levels = levels)
#df$FID2.type <- factor(df$FID2.type, levels = levels)
#df$Species[is.na(df$Species)] <- "None"
#df[is.na(df$FID1.type),'FID1.type'] <- "control"

#df[is.na(df$FID2.type),'FID2.type'] <- "control"


#head(df)

#table(df[df$FID1.type == 'control' & df$FID2.type == 'control',]$InfType)
#table(df[df$FID1.type == 'LR' & df$FID2.type == 'LR',]$InfType)
#write.csv(df[,c(2,4,1,5,3)],"FinalKing_JG/02.mergeQC/Final.type-match.usingAllSNP2ndDegree.csv",row.names = F)

#table(df[(df$FID1.type == 'control') | (df$FID2.type == 'control'),]$InfType)

#df[df$InfType == 'Dup/MZ',]
result_pair <-df
###########
head(real_pair)
colnames(real_pair)<-c("KCHIP_ID1","KCHIP_ID2","Rela_Pair")
head(sample.info)
df <- merge(real_pair,sample.info,by.y = "NewID",by.x = "KCHIP_ID1",all.x = T)
head(df)
colnames(df)[4] <- "KCHIP_ID1.type"
df <- merge(df,sample.info,by.y = "NewID",by.x = "KCHIP_ID2",all.x = T)
colnames(df)[5] <- "KCHIP_ID2.type"

df$KCHIP_ID1.type <- as.factor(df$KCHIP_ID1.type)
df$KCHIP_ID2.type <- as.factor(df$KCHIP_ID2.type)
real_pair <- df
########### 
head(result_pair)
str(result_pair)

head(result_pair)
head(real_pair)

result_pair$FID1 <- as.character(result_pair$FID1)
result_pair$FID2 <- as.character(result_pair$FID2)

head(result_pair)
head(real_pair)
tail(real_pair)

df <- merge(real_pair,result_pair,by.x = "KCHIP_ID1",by.y = "FID1",all.y = T)
#df <- merge(real_pair,result_pair,by.x = "KCHIP_ID1",by.y = "FID1",all.x = T)
head(df)
df1 <- merge(real_pair,result_pair,by.x = "KCHIP_ID2",by.y = "FID2",all.y = T)
#df1 <- merge(real_pair,result_pair,by.x = "KCHIP_ID2",by.y = "FID2",all.x = T)
df2 <- merge(df1,df,all = T)

########수정해야하는 사항
## NIH19KT6489 - NIH19KT6745 : spouse(5) -> FS(4)
table(df2$Rela_Pair)
df2[which(df2$KCHIP_ID1 == 'NIH19KT6489'),]$Rela_Pair <- "4"
str(df2)
table(df2$Rela_Pair)

df2[which(df2$KCHIP_ID1 == 'NIH19KT6489'),]
colnames(df2)
df2 <- as.data.frame(df2)
#########
head(df2)
ncol(df2)
table(df2$Rela_Pair)
#df2$Rela_Pair <- factor(df2$Rela_Pair,labels = c("PO","PO","MZ/twin","FS","Spouse","Relative","unrelated","Dual graft"))
df2$Rela_Pair <- factor(df2$Rela_Pair,labels = c("PO","Dual graft","PO","Dup/MZ","FS","Spouse","Relative","unrelated"))
#birthwt$smoke <- factor(birthwt$smoke,label = c("Non Smoker","Smoker"))
head(df)
head(df2)



df2 <-df2[,c("KCHIP_ID1","KCHIP_ID2","KCHIP_ID1.type","KCHIP_ID2.type","Rela_Pair","InfType","FID1","FID2","FID1.type","FID2.type")]
### type = type
#df <- df2[which(df2$FID1.type == df2$FID2.type),]
#head(df)
#write.csv(df,"04.Result/Final.kingResult.allsample.2ndDegree.allsnp.campare.clinical.and.king_same.type.csv",row.names = F,quote = F)
#write.csv(df2,"04.Result/Final.kingResult.allsample.2ndDegree.allsnp.campare.clinical.and.king.csv",row.names = F,quote = F)



########
#df <- df2[which(df2$FID1.type != df2$FID2.type),]
df <- df2
head(df)
df$Na <- 0
df[which(is.na(df$FID1) |is.na(df$FID2)|is.na(df$KCHIP_ID1)|is.na(df$KCHIP_ID2)),]$Na<-1
table(df$Na)

df$match <- 0
df[which(df$Rela_Pair == df$InfType),]$match <- 1
#df[-which(df$Rela_Pair == df$InfType),]$match <- 0
head(df)

#df[-which(is.na(df$FID1) |is.na(df$FID2)),]$Na <- 0
df[which(df$Na == 1),]$match <- 3

df[(df$KCHIP_ID1 == df$FID1 & df$KCHIP_ID2 == df$FID2) & df$match == 0,]$match <- 2
df[(df$KCHIP_ID1 == df$FID1 & df$KCHIP_ID2 == df$FID2) & df$match == 2  & df$Rela_Pair == 'Dual graft',]$match <- 6
df[(df$KCHIP_ID1 == df$FID1 & df$KCHIP_ID2 == df$FID2) & df$match == 2  & df$Rela_Pair == 'Relative',]$match <- 6
df[!(df$KCHIP_ID1 == df$FID1 & df$KCHIP_ID2 == df$FID2) & df$match == 1,]$match <- 5

df[(!is.na(df$KCHIP_ID1) & !is.na(df$KCHIP_ID2)) & df$Na == 1,]$match <- 4

#df[!(df$KCHIP_ID2 == df$FID1 & df$KCHIP_ID1 == df$FID2) & df$match == 1,]$match <- 4

#write.csv(df,"04.Result/Final.kingResult.allsample.2ndDegree.allsnp.campare.clinical.and.king.csv",row.names = F,quote = F)
write.csv(df,"04.Result/Final.kingResult.allsample.2ndDegree.filter.snp.campare.clinical.and.king.csv",row.names = F,quote = F)
write.csv(df,"04.Result/Final.kingResult.allsample.2ndDegree.filter.pruning50.5.0.1.snp.campare.clinical.and.king.csv",row.names = F,quote = F)
#df <- read.csv("04.Result/Final.kingResult.allsample.2ndDegree.allsnp.campare.clinical.and.king.csv")
#head(df)



#####summary
df <- read.csv("04.Result/Final.kingResult.allsample.2ndDegree.allsnp.campare.clinical.and.king.csv")
df <- read.csv("04.Result/Final.kingResult.allsample.2ndDegree.filter.snp.campare.clinical.and.king.csv")
df <- read.csv("04.Result/Final.kingResult.allsample.2ndDegree.filter.pruning50.5.0.1.snp.campare.clinical.and.king.csv")
table(df$match)
table(df$InfType)
table(king$InfType)


##############
#######MZ
table(real$Rela_Pair)

real$Rela_Pair <- factor(real$Rela_Pair,labels = c("PO","PO","Dup/MZ","FS","Spouse","Relative","unrelated","Deceased Heart beating","Non-heart beating","Dual graft"))
table(real$Rela_Pair)
nrow(real)
real[real$Rela_Pair == "Dual graft",]
table(result_pair$InfType)

result_pair[result_pair$InfType == "Dup/MZ",]

######summary
table(king$InfType)
sum(table(king$InfType))


head(as.data.frame(king$FID2))
a <- merge(data.frame(king$FID1),data.frame(king$FID2),by.x = "king.FID1",by.y = "king.FID2",all = T)
head(a)


real_pair$Rela_Pair <- factor(real_pair$Rela_Pair,labels = c("PO","PO","MZ/twin","FS","Spuse","Relative","unrelated","Dual graft"))
table(real_pair$Rela_Pair)/2
sum(table(real_pair$Rela_Pair))/2

b <- merge(data.frame(real_pair$KCHIP_ID1),data.frame(real_pair$KCHIP_ID2),by.x = "real_pair.KCHIP_ID1",by.y = "real_pair.KCHIP_ID2",all = T)

####check number of pair
head(real_pair)
head(sample.info)
a <- merge(real_pair,sample.info,by.x = "KCHIP_ID",by.y = "NewID")
a <- merge(real_pair,sample.info,by.x = "KCHIP_ID1",by.y = "NewID")
a

nrow(real_pair[!is.na(real_pair$Rela_Pair1),])/2
table(real_pair$Rela_Pair1)/2
real_pair$Rela_Pair1 <- factor(real_pair$Rela_Pair1,labels = c("PO","PO","MZ/twin","FS","Spouse","Relative","unrelated","Dual graft"))
table(real_pair$Rela_Pair1)/2
table(real_pair$Rela_Pair1)[c("PO","MZ/twin","FS","Relative")]/2
sum(table(real_pair$Rela_Pair1)[c("PO","MZ/twin","FS","Relative")]/2)

a <- real[which(real$Rela_Pair == 10),]
a
