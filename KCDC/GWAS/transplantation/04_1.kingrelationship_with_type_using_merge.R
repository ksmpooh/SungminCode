#######################20200422 merge king Result procossing
#####ctrl - ctrl, LD - LD 등 성향 정리
setwd("c:/Users/user/Desktop/KCDC/")
#king <- read.table("FinalKing_JG/02.mergeQC/01.allsnp/2ndDegree/JG.all.sample.2ndDegree.allsnp.related.kin0",header = T)
king <- read.table("FinalKing_JG/02.mergeQC/02.pruning/DefaultDegree/JG.merge.rmMHC.prunedsnp.defaultDegree.kin0",header = T)
sample.info <- read.table("transplantation/final_sample_info/last.sample.info.txt",header = T)
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
levels <- levels(df$FID1.type)

levels[length(levels) + 1] <- "control"

# refactor Species to include "None" as a factor level
# and replace NA with "None"
df$FID1.type <- factor(df$FID1.type, levels = levels)
df$FID2.type <- factor(df$FID2.type, levels = levels)
#df$Species[is.na(df$Species)] <- "None"
df[is.na(df$FID1.type),'FID1.type'] <- "control"

df[is.na(df$FID2.type),'FID2.type'] <- "control"


head(df)

table(df[df$FID1.type == 'control' & df$FID2.type == 'control',]$InfType)
table(df[df$FID1.type == 'LR' & df$FID2.type == 'LR',]$InfType)
#write.csv(df[,c(2,4,1,5,3)],"FinalKing_JG/02.mergeQC/Final.type-match.usingAllSNP2ndDegree.csv",row.names = F)

table(df[(df$FID1.type == 'control') | (df$FID2.type == 'control'),]$InfType)

df[df$InfType == 'Dup/MZ',]




