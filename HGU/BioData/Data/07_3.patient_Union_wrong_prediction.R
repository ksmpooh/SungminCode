#틀린것 찾는 venn

library(limma)
library(gplots)
library(dplyr)

wdir <- "D:/biodatalab/2018-1/match_wrong_prediction/"

ordata <- read.csv(paste0(wdir,"whole_result.csv"))
rownames(ordata)<-ordata$patient
####
for(i in 0:4){
  df <- read.csv(paste0(wdir,"te_result_",i,".csv"))
  rownames(df)<- df$patient
  colnames(df)<-c("index","patient","prediction")
  assign(paste0("ts_result_",i),df)
  

  df_2 <- ordata[rownames(df),]
  df_2 <- subset(df_2, select = "result")
  mdf <- cbind(df,df_2)
  sub <-df[(mdf$prediction != mdf$result),] 
  assign(paste0("ts_wrong_rownames_",i) , rownames(sub))
  
  df <- read.csv(paste0(wdir,"tr_result_",i,".csv"))
  rownames(df)<- df$patient
  colnames(df)<-c("index","patient","prediction")
  assign(paste0("tr_result_",i),df)
  
  df_2 <- ordata[rownames(df),]
  df_2 <- subset(df_2, select = "result")
  mdf <- cbind(df,df_2)
  sub <-df[(mdf$prediction != mdf$result),] 
  assign(paste0("tr_wrong_rownames_",i) , rownames(sub))
  
}


##ts
udf <- union(ts_wrong_rownames_0,ts_wrong_rownames_1)
udf <- union(udf,ts_wrong_rownames_2)
udf <- union(udf,ts_wrong_rownames_2)
udf <- union(udf,ts_wrong_rownames_2)

all <- Reduce(rbind,udf)
all <- as.data.frame(all)

all <- t(all)
all <-as.data.frame(all)
colnames(all)<-udf

all$index <- "all"

ts_0 <-as.data.frame(ts_wrong_rownames_0)
ts_0 <- t(ts_0)
ts_0 <- as.data.frame(ts_0)
colnames(ts_0)<-ts_wrong_rownames_0
ts_0$index <- "ts_0"

ts_1 <-as.data.frame(ts_wrong_rownames_1)
ts_1 <- t(ts_1)
ts_1 <- as.data.frame(ts_1)
colnames(ts_1)<-ts_wrong_rownames_1
ts_1$index <- "ts_1"

ts_2 <-as.data.frame(ts_wrong_rownames_2)
ts_2 <- t(ts_2)
ts_2 <- as.data.frame(ts_2)
colnames(ts_2)<-ts_wrong_rownames_2
ts_2$index <- "ts_2"

ts_3 <-as.data.frame(ts_wrong_rownames_3)
ts_3 <- t(ts_3)
ts_3 <- as.data.frame(ts_3)
colnames(ts_3)<-ts_wrong_rownames_3
ts_3$index <- "ts_3"

ts_4 <-as.data.frame(ts_wrong_rownames_4)
ts_4 <- t(ts_4)
ts_4 <- as.data.frame(ts_4)
colnames(ts_4)<-ts_wrong_rownames_4
ts_4$index <- "ts_4"

df <-data.frame()

df <- bind_rows(all,ts_0)
df <- bind_rows(df,ts_1)
df <- bind_rows(df,ts_2)
df <- bind_rows(df,ts_3)
df <- bind_rows(df,ts_4)

rownames(df)<-df$index
df<-subset(df,select = -index)

df_ <- t(df)
df_ <-as.data.frame(df_)
#colnames(df_)
#rownames(df_)
df_ <-subset(df_,select = -all)

id <- (df_!="")
id[is.na(id)]<-FALSE 
id.df<-as.data.frame(id)

vennCounts(id.df)
venn(id.df)



##tr



udf <- union(tr_wrong_rownames_0,tr_wrong_rownames_1)
udf <- union(udf,tr_wrong_rownames_2)
udf <- union(udf,tr_wrong_rownames_2)
udf <- union(udf,tr_wrong_rownames_2)


all <- Reduce(rbind,udf)
all <- as.data.frame(all)

all <- t(all)
all <-as.data.frame(all)
colnames(all)<-udf

all$index <- "all"



tr_0 <-as.data.frame(tr_wrong_rownames_0)
tr_0 <- t(tr_0)
tr_0 <- as.data.frame(tr_0)
colnames(tr_0)<-tr_wrong_rownames_0
tr_0$index <- "tr_0"

tr_1 <-as.data.frame(tr_wrong_rownames_1)
tr_1 <- t(tr_1)
tr_1 <- as.data.frame(tr_1)
colnames(tr_1)<-tr_wrong_rownames_1
tr_1$index <- "tr_1"

tr_2 <-as.data.frame(tr_wrong_rownames_2)
tr_2 <- t(tr_2)
tr_2 <- as.data.frame(tr_2)
colnames(tr_2)<-tr_wrong_rownames_2
tr_2$index <- "tr_2"

tr_3 <-as.data.frame(tr_wrong_rownames_3)
tr_3 <- t(tr_3)
tr_3 <- as.data.frame(tr_3)
colnames(tr_3)<-tr_wrong_rownames_3
tr_3$index <- "tr_3"

tr_4 <-as.data.frame(tr_wrong_rownames_4)
tr_4 <- t(tr_4)
tr_4 <- as.data.frame(tr_4)
colnames(tr_4)<-tr_wrong_rownames_4
tr_4$index <- "tr_4"

df <-data.frame()

df <- bind_rows(all,tr_0)
df <- bind_rows(df,tr_1)
df <- bind_rows(df,tr_2)
df <- bind_rows(df,tr_3)
df <- bind_rows(df,tr_4)

rownames(df)<-df$index
df<-subset(df,select = -index)

df_ <- t(df)
df_ <-as.data.frame(df_)
#colnames(df_)
#rownames(df_)
df_ <-subset(df_,select = -all)

id <- (df_!="")
id[is.na(id)]<-FALSE 
id.df<-as.data.frame(id)

vennCounts(id.df)
venn(id.df)



