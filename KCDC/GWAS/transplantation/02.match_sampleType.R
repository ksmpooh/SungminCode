##################smample information by type

setwd("c:/Users/user/Desktop/KCDC/transplantation/")

jg <- read.csv("sample_info/KchipJG_pairTable_20190908.csv")
head(jg)
rec <- subset(jg,select = c("OriID","NewID"))
don <- subset(jg,select = c("OriID.1","NewID.1"))
nrow(jg) == nrow(rec) + nrow(don)

LR <- rec[grep("^LR",rec$OriID),]
KR <- rec[grep("^KR",rec$OriID),]
KD <- don[grep("^KD",don$OriID.1),]
LD <- don[grep("^LD",don$OriID.1),]


dim(LR) + dim(KR) + dim(LD) + dim(KD)


LR$type <- "LR"
KR$type <- "KR"
LD$type <- "LD"
KD$type <- "KD"

head(LD)

colnames(LD) <- c("OriID","NewID","type")
colnames(KD) <- c("OriID","NewID","type")

ref <- rbind(LR,KR)
ref <- rbind(ref,LD)
ref <- rbind(ref,KD)

write.table(ref,"sample_info/sample.info.with.type.txt",col.names = T,row.names = F,quote = F,sep = '\t')

head(ref)

############################3#############
rmPCA <- read.table("1stQC/rmlist/rmPCA.txt")
rmLQ <- read.table("1stQC/rmlist/rmLQSamples.txt")
colnames(rmLQ) <- c("NewID","NewID1")
colnames(rmPCA) <- c("NewID","NewID1")

rmLQ <- merge(rmLQ,ref,by = "NewID")

dim(rmLQ)
table(rmLQ$type)
head(rmLQ)


rmPCA <- merge(rmPCA,ref,by = "NewID")
dim(rmPCA)
table(rmPCA$type)


df <- rbind(rmLQ,rmPCA)
df <- df[,c(1,3,4)]
head(df)

write.table(df,"1st.rmlist.with.type.txt",col.names = T,row.names = F,quote = F,sep = '\t')





######last sample type
setwd("c:/Users/user/Desktop/KCDC/transplantation/")
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
dim(df)

dim(df)[1] - dim(liver)[1]*2 - dim(kidney)[1] * 2

table(df$type)


ori <- read.table("summary.info.txt")
ori <- ori[order(ori$V1),]
dim(ori)
head(ori)
ori <- ori[,c("V1","V5")]
colnames(ori) <- c("NewID","sex")
head(df)
dim(df)
a <- merge(ori,df,by = "NewID")
head(a)
dim(a)
table(a$type)
cel <- read.table("cel_file_list.txt",header = T)
head(cel)
library(stringr)
cel$NewID <- str_split_fixed(str_split_fixed(cel$cel_files,"_",6)[,6],".CEL",2)[,1]
head(cel)
dim(cel)
write.table(cel,"sample_info/celfiles.and.NIHID.6772.txt",col.names = T,row.names = F,quote = F)

df <- merge(a,cel,by = "NewID")
head(df)
dim(df)
levels(df$type)

head(df)
dim(df)

write.table(df,"sample_info/last.sample.info.txt",col.names = T,row.names = F,quote = F)



########add etc
df <- read.table("sample_info/last.sample.info.txt",header = T)
df <- df[,1:ncol(df)-1]
#df$state <- "use"
#df$info <- "Normal"
rmLQ <- read.table("1stQC/rmlist/rmLQSamples.txt")
rmPCA <- read.table("1stQC/rmlist/rmPCA.txt")
head(rmLQ)
head(rmPCA)

rmLQ$info <- "miss-het"
rmLQ$state <- "remove"
rmPCA$info <- "PCA"
rmPCA$state <- "remove"
rmlist <- rbind(rmLQ,rmPCA)
dim(rmlist)
head(rmlist)
rmlist <- rmlist[,2:4]
colnames(rmlist)[1]<-"NewID"
head(rmlist)
#df[df$NewID == rmlist$NewID,'info'] <- rmlist$info

df <- merge(df,rmlist,all.x = TRUE,by = 'NewID')
dim(df)
head(df)
table(df$state)
table(df$info)

#a <- df[1:3,]
#a$state <- 'remove'

#df1<-df
#df1[df1$NewID == a$NewID,"state"] <- a$state
#head(df1)

head(df)
king <- read.table("sample_info/kinginfo.txt",header = T)
dim(king)
head(king)

pair <- read.csv("KchipJG_pairTable_20190908.csv",header = T)
head(pair)
dim(pair)
table(pair$Match)
pair <- subset(pair,pair$Match == 'MATCH')
pair <- subset(pair,select = c("Match","NewID","NewID.1"))
pair01 <- pair
colnames(pair01)[3]<-"pairID"
pair02 <- pair[,c(1,3)]
colnames(pair02)[2]<-"NewID"

pair02$pairID <- pair01$NewID
head(pair02)
head(pair01)

pair <-rbind(pair01,pair02)
head(pair)
dim(pair)
df <- merge(df,pair,by = "NewID",all.x = TRUE )
dim(df)
head(df)
df[is.na(df$info),'info'] <- "normal"
df[is.na(df$state),'state'] <- "use"