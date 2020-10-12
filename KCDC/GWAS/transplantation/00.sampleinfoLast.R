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

df <- read.table("c:/Users/user/Desktop/KCDC/transplantation/final_sample_info/last.sample.info.txt",header = T)
table(df$state)


##일반용역 정보 추가

#### HLa 일반용역 sample 정리

setwd("c:/Users/user/Desktop/KCDC/")
df1 <- read.csv("일반용역/2020년 장기이식 유전체정보 생산 대상자(DNALink 보유 검체 중)_20200623.csv")
df2 <- read.csv("transplantation/HLAtyping/20200828/HLAtyping.alle.gene.2digit.csv")
head(df1)
head(df2)

ref <- read.table("transplantation/final_sample_info/last.sample.info.txt",header = T)
head(ref)
rownames(ref)<-ref$NewID
ref$HLAtyping <- 0

ref[df1$NIHID,]$HLAtyping <- 2020
ref[df2$KID,]$HLAtyping <- 2019

table(ref$HLAtyping)
head(ref)
colnames(ref)[1] <- "KBA_ID"
write.table(ref,"transplantation/final_sample_info/last.sample.info_20200929.txt",col.names = T,row.names = F,quote = F)


out <- merge(df1,ref,by.x = "NIHID",by.y = "NewID",all.x = T)
out$HLA.ID <- "2020KDCA"
for (i in 1:265) {
  if (i < 10) {
    out[i,'HLA.ID'] <-paste0("2020KDCA00",i)  
  } else if ( i < 100){
    out[i,'HLA.ID'] <-paste0("2020KDCA0",i)
  } else{
    out[i,'HLA.ID'] <-paste0("2020KDCA",i)
  }
}
head(out)
write.csv(out,"일반용역/2020.HLAtyping.265sample.original.info.csv",row.names=F,quote=F)

#### 1차년도 아이디와 matching 및 수정
setwd("c:/Users/user/Desktop/KCDC/")
df <- read.csv("일반용역/2020년 장기이식 유전체정보 생산 대상자(DNALink 보유 검체 중)_20200623.csv")
head(df)
ref <- read.table("transplantation/final_sample_info/last.sample.info_20200929.txt",header = T)
head(ref)
df2 <- read.csv("일반용역/2019년도HLAtyping.265.sample.csv",header = T)
head(df2)

df <- merge(df,ref[,c(1,3)],by.x="NIHID",by.y="KBA_ID")
head(df)

pairtable <- read.csv("transplantation/final_sample_info/한국인칩_장기이식_pairTable.csv")
head(pairtable)
pairtable <- pairtable[,c(2,3,8,9)]
head(pairtable)

a <- merge(df2,pairtable,by.x="생산아이디",by.y="NewID.1")
b <- merge(df,pairtable,by.x="NIHID",by.y="NewID")
head(a)
head(b)
b$OriID.x == b$OriID.y
a$OLD_ID == a$OriID.1
a <- a[,c(1,2,3,5,6)]
colnames(a) <- c("KBA_ID1","oriID1","HLAID1","KBA_ID2","oriID2")
head(a)
a$HLAID2 <- a$HLAID1


library(stringr)
for (i in 1:265) {
  a[i,"HLAID2"] <- str_replace(a[i,"HLAID2"],'H','2020KDCA')
}

ori <- read.csv("일반용역/2020.HLAtyping.265sample.original.info.csv",header=T)
head(ori)
head(a)
head(ref)
out <- merge(a,ori[,c(1,2,3)],by.x = "KBA_ID2",by.y="NIHID")
head(out)
write.csv(out,"일반용역/HLAtyping.265pairtable.csv",row.names=F,quote=F)

out <-merge(ori[,c(1:ncol(ori)-1)],a[,c(4,6)],by.x="NIHID",by.y = "KBA_ID2")
head(out)
write.csv(out,"일반용역/2020.HLAtyping.265sample.original.info.csv",row.names=F,quote=F)
