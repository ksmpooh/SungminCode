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

colnames(LD) <- c("OriID","NewID","type")
colnames(KD) <- c("OriID","NewID","type")

ref <- rbind(LR,KR)
ref <- rbind(ref,LD)
ref <- rbind(ref,KD)

write.table(ref,"sample_info/sample.info.with.type.txt",col.names = T,row.names = F,quote = F,sep = '\t')


############################3#############
