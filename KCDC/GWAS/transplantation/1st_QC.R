#############
library(stringr)

setwd("C:/Users/user/Desktop/KCDC/transplantation/2nd/")
df <- read.table("JG.list.txt")
head(df)
#df$ID <- str_split_fixed(df$cel_files,"_",6)[,6]
df$ID <- str_split_fixed(str_split_fixed(df$V1,"_",6)[,6],".CEL",2)[,1]

rmID <- read.table("rm.1stQC.list.txt")
row.names(df) <- df$ID

out <-df[rmID$V1,]
head(out)
write.table(out,"2nd_CEL_file_list.txt",row.names = F,col.names = F,quote = F,sep = "\t")

################3
rmpca <- read.table("../2nd/rmPCA.txt")
head(rmpca)
colnames(rmpca)<-c("NewID","NewID2")

LQ <- read.table("../2nd/rmLQSamples.txt")
colnames(LQ)<-c("NewID","NewID2")

jg.type <- read.table("../intensity/JG_type.txt",header = T)
head(jg.type)

rmpca <- merge(rmpca,jg.type,by = 'NewID' )
levels(rmpca$type)
summary(rmpca$type)


###################333
LQ <- merge(LQ,jg.type,by= 'NewID')
summary(LQ$type)

