setwd("c:/Users/user/Desktop/KCDC/transplantation/king/")


#df <- read.table("../king/test.kin0",header = T)


##########DATA result processing  ##################################
df <- read.table("2ndDegree/01.all.sample/01.all.snp/01.related/JG.all.snp.2nd.related.kin0",header = T)
df1 <- read.table("2ndDegree/01.all.sample/02.pruned.snp/01.related/JG.pruned.snp.2nd.related.kin0",header = T)

a <- as.data.frame(table(df$InfType))
a <- data.frame(t(as.data.frame(table(df$InfType))))
colnames(a) <- as.data.frame(table(df$InfType))$Var1
a <- a[2,]
a$type <- "allsample.allsnp"
rownames(a) <- NULL
#a$UN <- 'NA'


b <- as.data.frame(table(df1$InfType))
b <- data.frame(t(as.data.frame(table(df1$InfType))))
colnames(b) <- as.data.frame(table(df1$InfType))$Var1
b <- b[2,]
b$type <- "allsample.prunedsnp"
rownames(b) <- NULL
c <- merge(a,b,all=T)
#c <- rbind(a,b)
c
df <- read.table("2ndDegree/02.notpredict/01.all.snp/01.related/JG.notpredict.related.2nd.kin0",header = T)
df1 <- read.table("2ndDegree/02.notpredict/02.pruned.snp/01.related/JG.not.predict.pruned.related.2nd.kin0",header = T)

a <- as.data.frame(table(df$InfType))
a <- data.frame(t(as.data.frame(table(df$InfType))))
colnames(a) <- as.data.frame(table(df$InfType))$Var1
a <- a[2,]
a$type <- "notpredict.allsnp"
rownames(a) <- NULL
#a$UN <- 'NA'
a
c <- merge(c,a,all = T)
c

b <- as.data.frame(table(df1$InfType))
b <- data.frame(t(as.data.frame(table(df1$InfType))))
colnames(b) <- as.data.frame(table(df1$InfType))$Var1
b <- b[2,]
b$type <- "notpredict.prunedsnp"
rownames(b) <- NULL
c <- merge(c,b,all = T)
c <- c[,c(4,6,3,2,1,5,7)]
c
###############################################################################################################
# 데이터 정리한 곳에 type 추가
#setwd("c:/Users/user/Desktop/KCDC/transplantation/king/")
setwd("c:/Users/user/Desktop/KCDC/FinalKing_JG/")
#df <- read.table("2ndDegree/01.all.sample/01.all.snp/01.related/JG.all.snp.2nd.related.kin0",header = T)
df <- read.table("02.mergeQC/01.allsnp/2ndDegree/JG.all.sample.2ndDegree.allsnp.related.kin0",header = T)
#df <- read.table("JG.2nd.merge.all.kin0",header = T)
head(df)
table(df$InfType)

df <- df[df$InfType == 'PO' | df$InfType == 'FS'|df$InfType =='Dup/MZ'|df$InfType =='2nd',c(1,3,ncol(df))]
#df <- df[df$InfType == 'PO' | df$InfType == 'FS'|df$InfType =='Dup/MZ',c(1,3,ncol(df))]
sample_info <- read.table("../transplantation/final_sample_info/last.sample.info.txt",header = T)
#sample_info <- read.table("../../final_sample_info/last.sample.info.txt",header = T)
head(sample_info)
sample_info <- sample_info[,c(1,2,4)]
head(df)

#a <- merge(sample_info,df,by.y = 'FID1',by.x = 'NewID')
#a <- merge(df,sample_info,by.x = 'FID1',by.y = 'NewID')
#colnames(a)[4:5] <- c("FID1.tubeID","FID1.type")
#head(a)

#b <- merge(df,sample_info,by.x = 'FID2',by.y = 'NewID')
#colnames(b)[4:5] <- c("FID2.tubeID","FID2.type")
#head(b)
#head(df)

#c <- merge(a,b,all = T)

c <- df
################20200401
a <- merge(df,sample_info,by.x = 'FID1',by.y = 'NewID',all.x = T)
colnames(a)[4:5] <- c("FID1.tubeID","FID1.type")
head(a)

b <- merge(df,sample_info,by.x = 'FID2',by.y = 'NewID',all.x = T)
colnames(b)[4:5] <- c("FID2.tubeID","FID2.type")

#a <- a[,c(1,2,5)]
c <- merge(a,b,all = T)


########################up 20200401

################################최초 중복 제거 위한 코드###########################3
## 다시 한 이유 : control 값이 사라짐
#c <- merge(b[,c(1,2,4,5),],c,by = 'FID2')
#d <- merge(b[,c(1,2,4,5),],d,by = 'FID2')
#d<-c

#head(c)


#library(dplyr)

#c <- distinct(c,FID1.x,FID2,.keep_all = TRUE) #중복된 값 제거
#head(d)
#c <- c[,c(1,2,3,4,6)]
#colnames(c)[2] <- 'FID1'
#head(c)

#c <- merge(c,a[,c(1,2,4,5)],by = 'FID1',all.x = TRUE)
#c <- distinct(c,FID2.x,FID1,.keep_all = TRUE)
################################last line#####최초 중복 제거 위한 코드###########################3


levels <- levels(c$FID1.type)
levels[length(levels) + 1] <- "control"

# refactor Species to include "None" as a factor level
# and replace NA with "None"
c$FID1.type <- factor(c$FID1.type, levels = levels)
c$FID2.type <- factor(c$FID2.type, levels = levels)
#df$Species[is.na(df$Species)] <- "None"
c[is.na(c$FID1.type),'FID1.type'] <- "control"

c[is.na(c$FID2.type),'FID2.type'] <- "control"
head(c)

out <- c[,c(1,5,4,2,7,6,3)]
#out <- c[,c(1,8,2,4,5,7,3)]
head(out)
#colnames(out) <- c('FID1','FID1.type','FID2','FID2.type','InfType','FID1.tubeID','FID2.tubeID')

out[out$FID2 == "NIH19KT0227",]

#out[grep("6715",out$FID1),]
out


#write.table(out,"Final/Final.kingResult.allsample.2ndDegree.allsnp.withType.txt",col.names = T,row.names = F,quote = F,sep = '\t')
#write.csv(out,"Final/Final.kingResult.allsample.2ndDegree.allsnp.withType.csv",row.names = F)
write.table(out,"03.FinalMatch/Final.kingResult.allsample.2ndDegree.allsnp.withType.txt",col.names = T,row.names = F,quote = F,sep = '\t')
write.csv(out,"03.FinalMatch/Final.kingResult.allsample.2ndDegree.allsnp.withType.csv",row.names = F)

dim(out)


