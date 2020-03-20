setwd("c:/Users/user/Desktop/KCDC/transplantation/2nd/all.king/")

#df <- read.table("../king/test.kin0",header = T)
df <- read.table("JG.2nd.merge.all.kin0",header = T)
head(df)
table(df$InfType)


df <- df[df$InfType == 'PO' | df$InfType == 'FS'|df$InfType =='Dup/MZ',c(1,3,ncol(df))]

sample_info <- read.table("../../final_sample_info/last.sample.info.txt",header = T)
head(sample_info)
sample_info <- sample_info[,c(1,2,4)]
head(df)

#a <- merge(sample_info,df,by.y = 'FID1',by.x = 'NewID')
a <- merge(df,sample_info,by.x = 'FID1',by.y = 'NewID')
colnames(a)[4:5] <- c("FID1.tubeID","FID1.type")
head(a)

b <- merge(df,sample_info,by.x = 'FID2',by.y = 'NewID')
colnames(b)[4:5] <- c("FID2.tubeID","FID2.type")
head(b)
head(df)

c <- df
c <- merge(b[,c(1,2,4,5),],c,by = 'FID2')
#d <- merge(b[,c(1,2,4,5),],d,by = 'FID2')

################################
################################
library(dplyr)

c <- distinct(c,FID1.x,FID2,.keep_all = TRUE) #중복된 값 제거
c <- c[,c(1,2,3,4,6)]
colnames(c)[2] <- 'FID1'
head(c)

c <- merge(c,a[,c(1,2,4,5)],by = 'FID1',all.x = TRUE)
c <- distinct(c,FID2.x,FID1,.keep_all = TRUE)

levels <- levels(c$FID1.type)
levels[length(levels) + 1] <- "control"

# refactor Species to include "None" as a factor level
# and replace NA with "None"
c$FID1.type <- factor(c$FID1.type, levels = levels)
#df$Species[is.na(df$Species)] <- "None"
c[is.na(c$FID1.type),'FID1.type'] <- "control"
head(c)

out <- c[,c(1,8,2,4,5,7,3)]
head(out)
colnames(out) <- c('FID1','FID1.type','FID2','FID2.type','InfType','FID1.tubeID','FID2.tubeID')

out[out$FID1 == "NIH19KT0016",]




write.table(out,"2ndkingResult.withType.txt",col.names = T,row.names = F,quote = F,sep = '\t')
dim(out)


