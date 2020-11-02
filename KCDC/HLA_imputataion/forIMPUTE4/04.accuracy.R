##################### compare result accuracy check

#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/HAN.ref/Result/")
#digit = "2"
#digit = "4"
#df <- read.csv("compare.IMPvsNGS.all.gene.2digit.csv",header = T)
#df <- read.csv("compare.IMPvsNGS.all.gene.4digit.csv",header = T)
#head(df)

#df<-df[!df$YSample == 'CDC015',]

check <- c("A.match","A.wrong","A.empty","B.match","B.wrong","B.empty","C.match","C.wrong","C.empty",
           "DRB1.match","DRB1.wrong","DRB1.empty","DPA1.match","DPA1.wrong","DPA1.empty",
           "DPB1.match","DPB1.wrong","DPB1.empty","DQA1.match","DQA1.wrong",
           "DQA1.empty","DQB1.match","DQB1.wrong","DQB1.empty")
length(check)
gene <- c("A","B","C","DRB1","DPA1","DPB1","DQA1","DQB1")
out.subset<-c("type","digit","A.match_SUM","A.wrong_SUM","A.empty_SUM","B.match_SUM","B.wrong_SUM","B.empty_SUM","C.match_SUM","C.wrong_SUM","C.empty_SUM",
              "DRB1.match_SUM","DRB1.wrong_SUM","DRB1.empty_SUM","DPA1.match_SUM","DPA1.wrong_SUM","DPA1.empty_SUM",
              "DPB1.match_SUM","DPB1.wrong_SUM","DPB1.empty_SUM","DQA1.match_SUM","DQA1.wrong_SUM",
              "DQA1.empty_SUM","DQB1.match_SUM","DQB1.wrong_SUM","DQB1.empty_SUM",
              "A_accuracy","B_accuracy","C_accuracy","DRB1_accuracy","DPA1_accuracy","DPB1_accuracy","DQA1_accuracy","DQB1_accuracy",
              "accuracy","accuracy(DQA1.out)")


accuracy.cal <- function(df,type,digit,check,gene,out.subset){
  ref <- matrix(nrow = 1,ncol = length(out.subset))
  ref <- as.data.frame(ref)
  colnames(ref)<-out.subset
  ref[,"type"] <- type
  ref[,"digit"] <- digit
  #print(ref)
  
  for (i in 1:length(check)){
    ref[1,paste0(check[i],"_SUM")] <- sum(df[,check[i]])
  }
  for (i in 1:length(gene)){
    ref[1,paste0(gene[i],"_accuracy")] <- ref[1,paste0(gene[i],".match_SUM")]/(ref[1,paste0(gene[i],".match_SUM")]+ref[1,paste0(gene[i],".wrong_SUM")])
  }
  ref[1,"accuracy"] <- sum(ref[1,c("A_accuracy","B_accuracy","C_accuracy","DRB1_accuracy","DPA1_accuracy","DPB1_accuracy","DQA1_accuracy","DQB1_accuracy")])/8
  ref[1,"accuracy(DQA1.out)"] <- sum(ref[1,c("A_accuracy","B_accuracy","C_accuracy","DRB1_accuracy","DPA1_accuracy","DPB1_accuracy","DQB1_accuracy")])/7
  return(ref)
  #print(ref)
}


ref = "Pan"
ref = "Han"
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/PAN.ref/Result/")
#setwd(paste0("c:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/",ref,".ref/Result/"))
setwd(paste0("c:/Users/user/Desktop/KCDC/HLAimputation/20201026/IMPUTE4/",ref,"/"))

digit = "2"

df <- read.csv(paste0("compare.IMPvsNGS.all.gene.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

a <- accuracy.cal(df,paste0(ref,".impute4.sample6574"),digit,check,gene,out.subset)

digit = "4"
df <- read.csv(paste0("compare.IMPvsNGS.all.gene.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

b <- accuracy.cal(df,paste0(ref,".impute4.sample6574"),digit,check,gene,out.subset)

#a <- accuracy.cal(df,"Han.impute4","2",check,gene,out.subset)
a
b
out <- rbind(a,b)
out1 <- rbind(a,b)

out <- rbind(out,out1)
out

###############Á¤¸®

#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/02.han/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/01.pan/")
ref = "Han"
ref = "Pan"
digit = "2"

df <- read.csv(paste0("compare.IMPvsNGS.all.gene.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]
a <- accuracy.cal(df,paste0(ref,".cookHLA.sample255"),digit,check,gene,out.subset)

digit = "4"
df <- read.csv(paste0("compare.IMPvsNGS.all.gene.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

b <- accuracy.cal(df,paste0(ref,".cookHLA.sample255"),digit,check,gene,out.subset)
out <- rbind(out,a)
out <- rbind(out,b)
out$type
out

###################
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20200731/")
ref = 'Han'
ref = 'Pan'
digit = "2"

df <- read.csv(paste0(ref,"/compare.IMPvsNGS.all.gene.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]
head(df)
a <- accuracy.cal(df,paste0(ref,".cookHLA.sample6574"),digit,check,gene,out.subset)

digit = "4"
df <- read.csv(paste0(ref,"/compare.IMPvsNGS.all.gene.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

b <- accuracy.cal(df,paste0(ref , ".cookHLA.sample6574"),digit,check,gene,out.subset)
head(a)
head(b)
#out <- read.csv("c:/Users/user/Desktop/KCDC/HLAimputation/cookHLAvsIMPUTE4.compare.Resul t.csv",header = T)
out <- rbind(out,a)
out <- rbind(out,b)
out$type
out$digit

write.csv(out,"c:/Users/user/Desktop/KCDC/HLAimputation/cookHLAvsIMPUTE4.compare.Result.csv",col.names = T,row.names = F,quote = F)


############20201102 A,B DRB1 

check <- c("A.match","A.wrong","A.empty","B.match","B.wrong","B.empty",
           "DRB1.match","DRB1.wrong","DRB1.empty")
length(check)
gene <- c("A","B","DRB1")
out.subset<-c("type","digit","A.match_SUM","A.wrong_SUM","A.empty_SUM","B.match_SUM","B.wrong_SUM","B.empty_SUM",
              "DRB1.match_SUM","DRB1.wrong_SUM","DRB1.empty_SUM",
              "A_accuracy","B_accuracy","DRB1_accuracy",
              "accuracy")



accuracy.cal <- function(df,type,digit,check,gene,out.subset){
  ref <- matrix(nrow = 1,ncol = length(out.subset))
  ref <- as.data.frame(ref)
  colnames(ref)<-out.subset
  ref[,"type"] <- type
  ref[,"digit"] <- digit
  #print(ref)
  
  for (i in 1:length(check)){
    ref[1,paste0(check[i],"_SUM")] <- sum(df[,check[i]])
  }
  for (i in 1:length(gene)){
    ref[1,paste0(gene[i],"_accuracy")] <- ref[1,paste0(gene[i],".match_SUM")]/(ref[1,paste0(gene[i],".match_SUM")]+ref[1,paste0(gene[i],".wrong_SUM")])
  }
  ref[1,"accuracy"] <- sum(ref[1,c("A_accuracy","B_accuracy","DRB1_accuracy")])/3
  #ref[1,"accuracy(DQA1.out)"] <- sum(ref[1,c("A_accuracy","B_accuracy","C_accuracy","DRB1_accuracy","DPA1_accuracy","DPB1_accuracy","DQB1_accuracy")])/7
  return(ref)
  #print(ref)
}

ref = "Pan"
ref = "Han"
setwd(paste0("c:/Users/user/Desktop/KCDC/HLAimputation/20201026/IMPUTE4/",ref,"/"))

digit = "2"

df <- read.csv(paste0("compare.IMPvsNGS.A.B.DRB1.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

a <- accuracy.cal(df,paste0(ref,".impute4"),digit,check,gene,out.subset)

digit = "4"
df <- read.csv(paste0("compare.IMPvsNGS.A.B.DRB1.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

b <- accuracy.cal(df,paste0(ref,".impute4"),digit,check,gene,out.subset)

#a <- accuracy.cal(df,"Han.impute4","2",check,gene,out.subset)
a
b
out <- rbind(a,b)
out1 <- rbind(a,b)

out <- rbind(out,out1)
out



ref = "Pan"
ref = "Han"
setwd(paste0("c:/Users/user/Desktop/KCDC/HLAimputation/20201026/cookHLA/",ref,"/"))

digit = "2"

df <- read.csv(paste0("compare.IMPvsNGS.A.B.DRB1.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

a <- accuracy.cal(df,paste0(ref,".cookHLA"),digit,check,gene,out.subset)

digit = "4"
df <- read.csv(paste0("compare.IMPvsNGS.A.B.DRB1.",digit,"digit.csv"),header = T)
df<-df[!df$YSample == 'CDC015',]

b <- accuracy.cal(df,paste0(ref,".cookHLA"),digit,check,gene,out.subset)


out <- rbind(out,a)
out <- rbind(out,b)

out



write.csv(out,"c:/Users/user/Desktop/KCDC/HLAimputation/20201026/cookHLAvsIMPUTE4.compare.Result.csv",col.names = T,row.names = F,quote = F)
``

