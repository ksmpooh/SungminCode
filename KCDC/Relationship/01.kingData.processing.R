####### 20201209 Final king

setwd("c:/Users/user/Desktop/KCDC/FinalKing_JG/")

ref <- read.table("../transplantation/final_sample_info/last.sample.info_20201207.txt",header = T)
head(ref)
real <- read.csv("04.newKing.6772/00.ref/KOTRY_KCHIP_Realpair.onlyPairInfo.20201209.csv")
head(real)
table(real$Rela_Pair)
#colnames(real)<-c("REAL.ID1","REAL.ID2","REAL.pair")
#real <- merge(real,ref[,c("KBA_ID","OriID","type")],by.x = "REAL.ID1",by.y="KBA_ID")
#colnames(real)[4:5]<-c("REAL.ID1.OriID","REAL.ID1.type")
#real <- merge(real,ref[,c("KBA_ID","OriID","type")],by.x = "REAL.ID2",by.y="KBA_ID")
#colnames(real)[6:7]<-c("REAL.ID2.OriID","REAL.ID2.type")
#real <- real[,c("REAL.ID1","REAL.ID1.OriID","REAL.ID1.type","REAL.ID2","REAL.ID2.OriID","REAL.ID2.type","REAL.pair")]

colnames(real)<-c("REAL.ID1","REAL.ID2","REAL.pair")
real <- merge(real,ref[,c("KBA_ID","type")],by.x = "REAL.ID1",by.y="KBA_ID")
colnames(real)[4]<-c("REAL.ID1.type")
real <- merge(real,ref[,c("KBA_ID","type")],by.x = "REAL.ID2",by.y="KBA_ID")
colnames(real)[5]<-c("REAL.ID2.type")
#real <- real[,c("REAL.ID1","REAL.ID1.type","REAL.ID2","REAL.ID2.type","REAL.pair")]
real <- real[,c("REAL.ID1","REAL.ID2","REAL.ID1.type","REAL.ID2.type","REAL.pair")]

## king결과에 oriID, type 추가
preprocessing_kingResult.with.type <- function(datain,refin){
  datain <- subset(datain,select = c("FID1","FID2","InfType"))
  colnames(datain)[1:2] <- c("Inf.ID1","Inf.ID2")
  datain <- merge(datain,refin[,c("KBA_ID","OriID","type")],by.x = "Inf.ID1",by.y="KBA_ID")
  colnames(datain)[4:5]<-c("Inf.ID1.OriID","Inf.ID1.type")
  datain <- merge(datain,refin[,c("KBA_ID","OriID","type")],by.x = "Inf.ID2",by.y="KBA_ID")
  colnames(datain)[6:7]<-c("Inf.ID2.OriID","Inf.ID2.type")
  return(datain[,c("Inf.ID1","Inf.ID1.OriID","Inf.ID1.type","Inf.ID2","Inf.ID2.OriID","Inf.ID2.type","InfType")])
}
preprocessing_kingResult.with.type <- function(datain,refin){
  datain <- subset(datain,select = c("FID1","FID2","InfType"))
  colnames(datain)[1:2] <- c("Inf.ID1","Inf.ID2")
  datain <- merge(datain,refin[,c("KBA_ID","type")],by.x = "Inf.ID1",by.y="KBA_ID")
  colnames(datain)[4]<-c("Inf.ID1.type")
  datain <- merge(datain,refin[,c("KBA_ID","type")],by.x = "Inf.ID2",by.y="KBA_ID")
  colnames(datain)[5]<-c("Inf.ID2.type")
  #return(datain[,c("Inf.ID1","Inf.ID1.type","Inf.ID2","Inf.ID2.type","InfType")])
  return(datain[,c("Inf.ID1","Inf.ID2","Inf.ID1.type","Inf.ID2.type","InfType")])
}



## pair match
# pair_match(kingdf,realdf,string="allsnp",integer=1)
# choice == 1 단순 페어 매치수
# choice == 2 이면 페어된 data frame
pair_match<-function(kingin,realin,concept,choice){
  colnames(kingin)[1:2]<- c("KCHIPID1","KCHIPID2")
  colnames(realin)[1:2]<- c("KCHIPID1","KCHIPID2")
  out <- merge(kingin,realin,by=c("KCHIPID1","KCHIPID2"))
  if (choice == 1) {
    a <- as.data.frame(table(out$REAL.pair == out$InfType))
    df<-data.frame(type = concept,pair_match = nrow(out),type_match =a[a$Var1 == 'TRUE',]$Freq)  
    
    b <- out[out$InfType == out$REAL.pair,]
    b <- as.data.frame(table(b$InfType))
    b <- as.data.frame(t(b))
    colnames(b) <- b["Var1",]
    b <- b['Freq',]
    rownames(b)<-1
    df <- merge(df,b)
    
    return(df)
  }else{
    return(out[,c(1,2,5,8)])  
  }
}

## king result
king.type.table<-function(kingin,concept){
  a <-as.data.frame(table(kingin$InfType))
  a <- t(a)    
  a <- as.data.frame(a)
  rownames(a)<-c(1,2)
  colnames(a)<-a[1,]
  a <- a[2,]
  a$type <- concept
  return(a)
}

df <- read.table("04.newKing.6772/01.allsnp/02.related.2ndDegree/related.2ndDegree.Result.kin0",header = T)
head(df)


df <- preprocessing_kingResult.with.type(df,ref)
head(df)
head(real)
a <- pair_match(df,real,"allsnp",2)
a <- a[a$InfType != a$REAL.pair,]
a

table(real$REAL.pair)/2
#######king result data table
out <- data.frame()
df <- read.table("04.newKing.6772/01.allsnp/02.related.2ndDegree/related.2ndDegree.Result.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- king.type.table(df,"allsnp")
df
out <-merge(out,df,all = T)
out <-merge(out,df,all = T)
out
df <- read.table("04.newKing.6772/02.filter/02.related.2ndDegree/filter.2ndDegree.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- king.type.table(df,"filter")
out <-merge(out,df,all = T)
out

df <- read.table("04.newKing.6772/03.pruning/02.related.2ndDegree/JG.pruning50.5.0.1.2ndDegree.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- king.type.table(df,"pruning")
out <-merge(out,df,all = T)

df <- read.table("04.newKing.6772/04.chr/JG.SNPolihser.chr1.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- king.type.table(df,"chr1")
out <-merge(out,df,all = T)

df <- read.table("04.newKing.6772/04.chr/JG.SNPolihser.fil.chr1.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- king.type.table(df,"filter.chr1")
out <-merge(out,df,all = T)

df <- read.table("04.newKing.6772/04.chr/JG.SNPolihser.fil.pruning.chr1.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- king.type.table(df,"pruning.chr1")
out <-merge(out,df,all = T)

out <- out[,c("type","Dup/MZ","PO","FS","2nd","3rd","4th","UN")]
out <- out[c(1,2,3,6,5,4),]
out


####### all type matching table 

out <- data.frame()
df <- read.table("04.newKing.6772/01.allsnp/02.related.2ndDegree/related.2ndDegree.Result.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- pair_match(df,real,"allsnp",1)
out<-rbind(out,df)
out

df <- read.table("04.newKing.6772/02.filter/02.related.2ndDegree/filter.2ndDegree.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- pair_match(df,real,"filter",1)
out<-rbind(out,df)
out
df <- read.table("04.newKing.6772/03.pruning/02.related.2ndDegree/JG.pruning50.5.0.1.2ndDegree.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- pair_match(df,real,"pruning",1)
out<-rbind(out,df)


df <- read.table("04.newKing.6772/04.chr/JG.SNPolihser.chr1.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- pair_match(df,real,"chr1",1)
out<-rbind(out,df)

df <- read.table("04.newKing.6772/04.chr/JG.SNPolihser.fil.chr1.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- pair_match(df,real,"filter.chr1",1)
out<-rbind(out,df)

df <- read.table("04.newKing.6772/04.chr/JG.SNPolihser.fil.pruning.chr1.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- pair_match(df,real,"pruning.chr1",1)
out<-rbind(out,df)
out
out <- out[,c(1,2,3,4,6,5)]

##test
df <- read.table("04.newKing.6772/01.allsnp/02.related.2ndDegree/related.2ndDegree.Result.kin0",header = T)
df <- preprocessing_kingResult.with.type(df,ref)
df <- pair_match(df,real,"allsnp",2)
head(df)
table(df[df$InfType != df$REAL.pair,]$REAL.pair)
table(df[df$InfType != df$REAL.pair,]$InfType)
head(real)
head(df)
a <- real[!(real$REAL.ID1 %in% df$KCHIPID1),]
table(a$REAL.pair)
table(real$REAL.pair)/2
table(df$InfType)

a <- df[df$InfType != df$REAL.pair,]
head(a)
table(a$InfType)
table(a$REAL.pair)
