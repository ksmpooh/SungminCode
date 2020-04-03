setwd("c:/Users/user/Desktop/KCDC/transplantation/")


library(dplyr)
#df <- read.table("king/data.processing/2ndkingResult.2ndDegree.allsnp.withType.txt",header = T)
#df <- read.table("king/Final/Final.kingResult.allsample.2ndDegree.allsnp.withType.txt",header = T)
#df <- read.table("")
#head(df)
#df[grep("6715",df$FID2),]
notpredict.allsnp <- read.table("king/2ndDegree/03.Final/Fianl.JG.not.predected.allsnp.2ndDegree.kin0",header = T)
head(notpredict.allsnp)
notpredict.allsnp[grep("6715",notpredict.allsnp$FID2),]
notpredict.allsnp <- notpredict.allsnp[notpredict.allsnp$InfType == 'PO'|notpredict.allsnp$InfType == 'FS'|notpredict.allsnp$InfType == '2nd',c(1,3,14)]
notpredict.allsnp


df[df$FID1 %in% notpredict.allsnp$FID1,c(1,3,5)]
#notpredict.allsnp


ori <- read.table("2nd/all.king/2ndkingResult.withType.txt",header = T)
#ori <- read.table("king/Final/Final.kingResult.allsample.2ndDegree.allsnp.withType.txt",header = T)
head(ori)
table(ori$InfType)


head(notpredict.list)
notpredict.list <- read.csv("king/compare.not.predict/notpredict.po.list.csv",header = T)
colnames(notpredict.list)[1]<- "KCHIP_ID"

head(notpredict.list)
notpredict.list <- rbind(notpredict.list,read.csv("king/compare.not.predict/notpredict.FS.list.csv",header = T))
notpredict.list <- rbind(notpredict.list,read.csv("king/compare.not.predict/notpredict.3rd.list.csv",header = T))
head(notpredict.list)
table(notpredict.list$Rela_Pair1)
notpredict.list
notpredict.allsnp
head(notpredict.allsnp)

#a<-merge(notpredict.allsnp,notpredict.list,by.y = "KCHIP_ID",by.x = "FID1",all = T)
#b<-merge(notpredict.allsnp,notpredict.list,by.y = "KCHIP_ID",by.x = "FID2")
a<-merge(notpredict.list,notpredict.allsnp,by.x = "KCHIP_ID",by.y = "FID1")
head(a)
b<-merge(notpredict.list,notpredict.allsnp,by.x = "KCHIP_ID",by.y = "FID2")
b
#colnames(b)[7] <- "FID2"
c <- merge(a,b,all = T)
head(c)



str(c)
c$Match <- 'NA'
#c$match <-lapply(c,function(x) ifelse ((x[,"FID1"] == x[,"KCHIP_ID1"]) | (x[,"FID2"] == x[,"KCHIP_ID"]),'Yes','NO'))
c$FID1 <- as.character(c$FID1)
c$FID2 <- as.character(c$FID2)
c$KCHIP_ID1 <- as.character(c$KCHIP_ID1)
c[is.na(c$FID2),]$FID2 <- 'NA'
c[is.na(c$FID1),]$FID1 <- 'NA'
#c <- c %>% mutate(FID1 = replace(FID1,is.na(FID1),0))
head(c)
c$match <-lapply((c$FID1 == c$KCHIP_ID1) | (c$FID2 == c$KCHIP_ID1),ifelse,"match","no")
table(c$match)


c[c$KCHIP_ID %in% notpredict.allsnp$FID1,]
notpredict.allsnp$FID1 %in% c$KCHIP_ID
notpredict.allsnp$FID1 %in% c$KCHIP_ID1

#write.csv(notpredict.allsnp[!(notpredict.allsnp$FID1 %in% c$KCHIP_ID) & !(notpredict.allsnp$FID1 %in% c$KCHIP_ID1),],"king/compare.not.predict/test.csv",row.names = F)

c[c$KCHIP_ID == "NIH19KT0010",]


c$match <- as.character(c$match)
head(c)
write.csv(c,"king/compare.not.predict/compare.notpredictlist.by.all.sample.king.with.king.Result,by.Notpredictlist.csv",row.names = F)



rmkinglist <- c("NIH19KT0595","NIH19KT0884","NIH19KT6929","NIH19KT7288","NIH19KT3613","NIH19KT4396")
d <- notpredict.list[notpredict.list$KCHIP_ID %in% rmkinglist,]
d$InfType = "rmkinglist(unrealtedSample).pair"
e <- notpredict.list[notpredict.list$KCHIP_ID1 %in% rmkinglist,]
e$InfType = "rmkinglist(unrealtedSample).pair"
e
head(c)
e
d <- merge(notpredict.list,e)
d

c<- merge(c,d,all = T)
c
#write.csv(c,"king/compare.not.predict/compare.notpredictlist.by.all.sample.king.with.king.Result,by.Notpredictlist.csv",row.names = F)

########################rmlist
rmlist <- read.table("final_sample_info/last.sample.info.txt",header = T)
head(rmlist)
rmlist <- subset(rmlist,rmlist$state == "remove",select = c("NewID","info","state"))


head(c)
d <- merge(c,rmlist,by.x = "KCHIP_ID",by.y = 'NewID',all.x = T)
head(d)
d$state


######################################################################20200402
#### FINAL #####
setwd("c:/Users/user/Desktop/KCDC/transplantation/")
library(dplyr)

#df <- read.table("king/data.processing/2ndkingResult.2ndDegree.allsnp.withType.txt",header = T)
#df <- read.table("king/Final/Final.kingResult.allsample.2ndDegree.allsnp.withType.txt",header = T)

notpredict.allsnp <- read.table("king/2ndDegree/03.Final/Fianl.JG.not.predected.allsnp.2ndDegree.kin0",header = T)
head(notpredict.allsnp)
notpredict.allsnp$InfType
#notpredict.allsnp <- notpredict.allsnp[notpredict.allsnp$InfType == 'PO'|notpredict.allsnp$InfType == 'FS'|notpredict.allsnp$InfType == '2nd',c(1,3,14)]
notpredict.allsnp <- notpredict.allsnp[,c(1,3,14)]
#notpredict.allsnp


#ori <- read.table("2nd/all.king/2ndkingResult.withType.txt",header = T)

notpredict.list <- read.csv("king/Final/Final.not.predicted.list.by.allsample.allsnp.king.csv",header = T)
notpredict.list <- subset(notpredict.list,select = c("KCHIP_ID","OriID","Rela_Pair","KCHIP_ID1","OriID1","Rela_Pair1"))
table(notpredict.list$Rela_Pair)
head(notpredict.list)
notpredict.list <- notpredict.list[(notpredict.list$Rela_Pair == 1 | notpredict.list$Rela_Pair == 2 | notpredict.list$Rela_Pair == 4 |notpredict.list$Rela_Pair == 6),]


head(notpredict.allsnp)

a<-merge(notpredict.list,notpredict.allsnp,by.x = "KCHIP_ID",by.y = "FID1")
head(a)
b<-merge(notpredict.list,notpredict.allsnp,by.x = "KCHIP_ID",by.y = "FID2")
b
#colnames(b)[7] <- "FID2"
c <- merge(a,b,all = T)
head(c)



str(c)
c$Match <- 'NA'
#c$match <-lapply(c,function(x) ifelse ((x[,"FID1"] == x[,"KCHIP_ID1"]) | (x[,"FID2"] == x[,"KCHIP_ID"]),'Yes','NO'))
c$FID1 <- as.character(c$FID1)
c$FID2 <- as.character(c$FID2)
c$KCHIP_ID1 <- as.character(c$KCHIP_ID1)
c[is.na(c$FID2),]$FID2 <- 'NA'
c[is.na(c$FID1),]$FID1 <- 'NA'
#c <- c %>% mutate(FID1 = replace(FID1,is.na(FID1),0))
head(c)
c$Match <-lapply((c$FID1 == c$KCHIP_ID1) | (c$FID2 == c$KCHIP_ID1),ifelse,"match","no")
table(c$Match)


c[c$KCHIP_ID %in% notpredict.allsnp$FID1,]
notpredict.allsnp$FID1 %in% c$KCHIP_ID
notpredict.allsnp$FID1 %in% c$KCHIP_ID1

#write.csv(notpredict.allsnp[!(notpredict.allsnp$FID1 %in% c$KCHIP_ID) & !(notpredict.allsnp$FID1 %in% c$KCHIP_ID1),],"king/compare.not.predict/test.csv",row.names = F)

c[c$KCHIP_ID == "NIH19KT0010",]


c$Match <- as.character(c$Match)
head(c)

#write.csv(c,"king/Final/Final.compare.notpredictlist.by.all.sample.king.with.king.Result,by.Notpredictlist.csv",row.names = F)



rmkinglist <- c("NIH19KT0595","NIH19KT0884","NIH19KT6929","NIH19KT7288","NIH19KT3613","NIH19KT4396")
d <- notpredict.list[notpredict.list$KCHIP_ID %in% rmkinglist,]
d$InfType = "rmkinglist(unrealtedSample).pair"
e <- notpredict.list[notpredict.list$KCHIP_ID1 %in% rmkinglist,]
e$InfType = "rmkinglist(unrealtedSample).pair"
e
head(c)
e
d <- merge(notpredict.list,e)
d

c<- merge(c,d,all = T)
c


write.csv(c,"king/Final/Final.compare.notpredictlist.by.all.sample.king.with.king.Result,by.Notpredictlist.csv",row.names = F)
