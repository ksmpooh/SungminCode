setwd("c:/Users/user/Desktop/KCDC/transplantation/")
cel <- read.table("final_sample_info/cel_file_list.txt",header = T)

library(stringr)
cel$NewID <- str_split_fixed(str_split_fixed(cel$cel_files,"_",6)[,6],".CEL",2)[,1]
write.table(cel,"final_sample_info/celfilelist.withNEw.txt",col.names = T,row.names = F,quote = F)

#####수정후....
cel <- read.table("final_sample_info/celfilelist.withNEw.txt",header = T)
ori <- read.csv("final_sample_info/한국인칩_장기이식_pairTable.csv",header = T)
head(ori)
donor <- ori[,8:10]
head(donor)
tail(donor)
donor <- donor[1:3034,]
rec <- ori[,2:4]
colnames(donor) <- colnames(rec)

ori <- rbind(donor,rec)
ori <- ori[,c(2,3)]
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
dim(df)

df <-merge(ori,df,by='NewID')
df <-merge(df,cel,by='NewID')

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





df[is.na(df$info),'info'] <- "normal"
df[is.na(df$state),'state'] <- "use"

write.table(df,"final_sample_info/last.sample.info.txt",col.names = T,row.names = F,quote = F)

df <- read.table("final_sample_info/last.sample.info.txt",header = T)
out <- df[grep("R",df$type),]
dim(out)
head(out)
write.table(out,"final_sample_info/donor.info.txt",col.names = T,row.names = F,quote = F)


#########################################
##20200812 추가
setwd("c:/Users/user/Desktop/KCDC/transplantation/final_sample_info/")
df<-read.table("last.sample.info.txt",header = T)
head(df)
table(df$info)
table(df$state)
df[which(df$info == 'PCA'),]$state <- '1st.remove'
df[which(df$info == 'miss-het'),]$state <- '1st.remove'

secondPCA <- read.table("../sampleQC/2nd.rmPCA.txt")
head(secondPCA)
head(df)
rownames(df) <- df$NewID
df[secondPCA$V1,]$info <- "PCA"
df[secondPCA$V1,]$state <- "2nd.remove"
table(df$state)
king <- read.table("../sampleQC/notRelatedId.byking.txt")
head(king)
df[king$V1,]$info <- "king"
df[king$V1,]$state <- "notRelated(MZ)"

table(df$state)
table(df$info)
head(df)
write.table(df,"last.sample.info.txt",col.names = T,row.names = F,quote = F,sep = "\t")
