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
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
dim(df)

dim(df)[1] - dim(liver)[1]*2 - dim(kidney)[1] * 2

table(df$type)


ori <- read.table("summary.info.txt")
dim(ori)
head(ori)
ori <- ori[,c("V1","V5")]
colnames(ori) <- c("NewID","sex")

a <- merge(ori,df,by = "NewID")
head(a)
dim(a)
table(a$type)
cel <- read.table("cel_file_list.txt",header = T)
head(cel)
library(stringr)
cel$ID <- str_split_fixed(str_split_fixed(cel$cel_files,"_",6)[,6],".CEL",2)[,1]
head(cel)

df <- merge(a,cel,by.x = 'NewID',by.y = 'ID')
head(df)
dim(df)
