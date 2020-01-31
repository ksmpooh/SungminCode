setwd("c:/Users/user/Desktop/KCDC/transplantation/intensity/")

library(stringr)
df <- read.table("test.txt",header = T)
dim(df)
head(df)
df$signal_sum_A_G <- df$axiom_signal_contrast_A_signal_mean + df$axiom_signal_contrast_G_signal_mean
df <- subset(df,select = c('cel_files','signal_sum_A_G'))
df$NewID <- str_split_fixed(df$cel_files,"_",6)[,6]

jg <- read.csv("../KchipJG_pairTable_20190908.csv",header = T)
head(jg)
dim(jg)
tail(jg)
rec <- subset(jg,select = c("OriID","NewID"))
don <- subset(jg,select = c("OriID.1","NewID.1"))
LR <- rec[grep("^LR",rec$OriID),]
KR <- rec[grep("^KR",rec$OriID),]
KD <- don[grep("^KD",don$OriID.1),]
LD <- don[grep("^LD",don$OriID.1),]

dim(LR)
dim(KR)
dim(LD)
dim(KD)
tail(LR)
tail(KR)
tail(KD)
tail(LD)
dim(LR) + dim(KR) + dim(LD) + dim(KD)


LR$type <- "LR"
KR$type <- "KR"
LD$type <- "LD"
KD$type <- "KD"

colnames(LD) <- c("OriID","NewID","type")
colnames(KD) <- c("OriID","NewID","type")
head(LR)
head(LD)
ref <- rbind(LR,KR)
ref <- rbind(ref,LD)
ref <- rbind(ref,KD)
head(ref)
dim(ref)
df <- merge(df,ref,by ="NewID" )
df1 <- merge(ref,df,by = "NewID")
dim(df1)
#df2 <- merge(df,ref,by ="NewID" )
dim(df2)
head(df)
dim(df)

boxplot(df$signal_sum_A_G~df$type,ylab = "signal_sum_A+G",xlab = "type"
        , main = "Transplantation Intensity")

###########################################
write.table(df,"JG_type.txt",col.names = T,row.names = F,quote = F,sep = '\t')

#################_2È®ÀÎ

tmp <- df[df$NewID != df2$NewID,]
head(tmp)
dim(tmp)
tmp <- setdiff(df$NewID,df2$NewID)
tmp
dim(df)
dim(df2)
df["NIH19KT4031_1","NewID"]
##############################3
jg <- read.csv("../KchipJG_pairTable_20190908.csv",header = T)
head(jg)
rec <- subset(jg,select = c("OriID","NewID"))
don <- subset(jg,select = c("OriID.1","NewID.1"))
LR <- rec[grep("^LR",rec$OriID),]
KR <- rec[grep("^KR",rec$OriID),]
KD <- don[grep("^KD",don$OriID.1),]
LD <- don[grep("^LD",don$OriID.1),]



