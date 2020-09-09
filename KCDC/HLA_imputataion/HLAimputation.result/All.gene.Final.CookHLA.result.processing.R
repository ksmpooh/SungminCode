### FINAL cookHLA result processing and compare with NGS

#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result3/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/UsingHan/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20200731/Han/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20200731/Pan/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/01.pan/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/02.han/")
#df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
head(df)
colnames(df)
#grep("*DRB1*",colnames(df))


library(stringr)

A <- df[,c(1,2,grep("HLA_A_*",colnames(df)))]
B <- df[,c(1,2,grep("HLA_B_*",colnames(df)))]
C <- df[,c(1,2,grep("HLA_C_*",colnames(df)))]
DRB1 <- df[,c(1,2,grep("*DRB1*",colnames(df)))]
#DRB3 <- df[,c(1,2,grep("*DRB3*",colnames(df)))]

DPB1 <- df[,c(1,2,grep("*DPB1*",colnames(df)))]
DPA1 <- df[,c(1,2,grep("*DPA1*",colnames(df)))]

DQB1 <- df[,c(1,2,grep("*DQB1*",colnames(df)))]
DQA1 <- df[,c(1,2,grep("*DQA1*",colnames(df)))]
head(DQA1)

hla.subset <- function(df,n){
  if( (n ==2) || (n == 4)){
    i = paste(n,"(nvalue) is OK")
    print(i)
    temp <- df[,1:2]
    #    print(temp)
    if(n == 2){
      for (i in 3:ncol(df)){
        #print(as.integer(str_split_fixed(colnames(df)[i],"_",4)[3]))
        if(as.integer(str_split_fixed(colnames(df)[i],"_",4)[3]) < 100){
          b <- df[,c(1,2,i)]
          temp <- merge(temp,b)
        }
      }
    }else{
      for (i in 3:ncol(df)){
        if(as.integer(str_split_fixed(colnames(df)[i],"_",4)[3]) >= 100){
          b <- df[,c(1,2,i)]
          temp <- merge(temp,b)
        }
      }
    }
    return(temp)
  }else{
    a = paste(n,"(nvalue) is wrong, only for 2,4")
    print(a)
    return(0)
  }
}
#DRB_td <- hla.subset(DRB,3)
hla.find<-function(df,concept){
  colcount <- ncol(df)
  print(colcount)
  rowcount <- nrow(df)
  print(rowcount)
  for(i in (1:rowcount)){
    n = 0
    for(j in (3:colcount)){
      if(df[i,j] == 2){
        if(n == 0){
          df[i,paste("IMP_",concept,'.1',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          df[i,paste("IMP_",concept,'.2',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 2  
        }else if(n == 1){
          df[i,paste("IMP_",concept,'.2',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          df[i,paste("IMP_",concept,'.3',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 2  
        }
      }else if(df[i,j] == 1){
        if(n == 0){
          df[i,paste("IMP_",concept,'.1',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 1
        }else if(n == 1){
          df[i,paste("IMP_",concept,'.2',sep = "")] <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
          n = n + 1
        }else{
          df[i,paste("IMP_",concept,'.3',sep = "")]  <- as.integer(str_split_fixed(colnames(df)[j],"_",4)[3])
        }
      }
    }
  }
  print(head(df))
  return(df)
}

#hla.subset(df,digit)
#hla.find(df_digit,"concept")



###################### 2digit
DRB1_td <- hla.subset(DRB1,2)
head(DRB1_td)
DRB1_td <- hla.find(DRB1_td,"DRB1")
table(DRB1_td$DRB3)
str(DRB1_td)

DQA1_td <- hla.subset(DQA1,2)
DQA1_td <- hla.find(DQA1_td,"DQA1")
DQB1_td <- hla.subset(DQB1,2)
DQB1_td <- hla.find(DQB1_td,"DQB1")


DPA1_td <- hla.subset(DPA1,2)
DPA1_td <- hla.find(DPA1_td,"DPA1")
DPB1_td <- hla.subset(DPB1,2)
DPB1_td <- hla.find(DPB1_td,"DPB1")


grep("*N_P",colnames(A))
colnames(A)[6]<-"HLA_A_0122_P"
colnames(A)[19]<-"HLA_A_0253_P"

A_td <- hla.subset(A,2)
A_td <- hla.find(A_td,"A")

B_td <- hla.subset(B,2)
B_td <- hla.find(B_td,"B")

C_td <- hla.subset(C,2)
C_td <- hla.find(C_td,"C")


table(DQA1_td$IMP_DQA1.3)
table(DQB1_td$IMP_DQB1.3)
table(DPA1_td$IMP_DPA1.3)
table(DPB1_td$IMP_DPB1.3)
table(DRB1_td$IMP_DRB1.3)
table(C_td$IMP_C.3)
table(B_td$IMP_B.3)
table(A_td$IMP_A.3)
head(A_td)

##HAN 2digit
DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID
##PAN 2digit
DPA1_td[!is.na(DPA1_td$IMP_DPA1.3),]$IID

out <- merge(A_td[,c('IID','IMP_A.1','IMP_A.2')],B_td[,c('IID','IMP_B.1','IMP_B.2')],by = 'IID')
out <- merge(out,C_td[,c('IID','IMP_C.1','IMP_C.2')],by = 'IID')
out <- merge(out,DRB1_td[,c('IID','IMP_DRB1.1','IMP_DRB1.2')],by = 'IID')
out <- merge(out,DPA1_td[,c('IID','IMP_DPA1.1','IMP_DPA1.2')],by = 'IID')
out <- merge(out,DPB1_td[,c('IID','IMP_DPB1.1','IMP_DPB1.2')],by = 'IID')
out <- merge(out,DQA1_td[,c('IID','IMP_DQA1.1','IMP_DQA1.2')],by = 'IID')
out <- merge(out,DQB1_td[,c('IID','IMP_DQB1.1','IMP_DQB1.2')],by = 'IID')


write.csv(out,"HLAimputation.all.gene.2digit.Result.csv",row.names = F,quote = F)

head(out)
#####merge NGS
out <- read.csv("HLAimputation.all.gene.2digit.Result.csv")
ngs <- read.csv("../../../transplantation/HLAtyping/20200828/HLAtyping.alle.gene.2digit.csv")
head(ngs)
ncol(ngs)
out <- merge(out,ngs,by.x = "IID",by.y = "KID")
ncol(out)
head(out)
out <-out[,c(1,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]
head(out)

#row.names(out) <- out$IID
#out[DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID,]




write.csv(out,"MERGE.impResult.hlatyping.all.gene.2digit.csv",row.names = F,quote = F)




#cmp.result.2digit <- read.csv("compare.IMPvsNGS.all.gene.2digit.csv")

#rownames(cmp.result.2digit) <- cmp.result.2digit$IID
#nrow(cmp.result.2digit[DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID,c("DQA1.match","DQA1.wrong","DQA1.empty")])






####################4 digit
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/01.pan/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/02.han/")
#df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
#df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
#head(df)
#colnames(df)
#grep("*DRB1*",colnames(df))


library(stringr)










DRB1_fd <- hla.subset(DRB1,4)
head(DRB1_fd)
DRB1_fd <- hla.find(DRB1_fd,"DRB1")
table(DRB1_fd$DRB3)


DQA1_fd <- hla.subset(DQA1,4)
DQA1_fd <- hla.find(DQA1_fd,"DQA1")
DQB1_fd <- hla.subset(DQB1,4)
DQB1_fd <- hla.find(DQB1_fd,"DQB1")


DPA1_fd <- hla.subset(DPA1,4)
DPA1_fd <- hla.find(DPA1_fd,"DPA1")
DPB1_fd <- hla.subset(DPB1,4)
DPB1_fd <- hla.find(DPB1_fd,"DPB1")


grep("*N_P",colnames(A))
colnames(A)[6]<-"HLA_A_0122_P"
colnames(A)[19]<-"HLA_A_0253_P"
A_fd <- hla.subset(A,4)
A_fd <- hla.find(A_fd,"A")

B_fd <- hla.subset(B,4)
B_fd <- hla.find(B_fd,"B")

C_fd <- hla.subset(C,4)
C_fd <- hla.find(C_fd,"C")


table(DQA1_fd$IMP_DQA1.3)
table(DQB1_fd$IMP_DQB1.3)
table(DPA1_fd$IMP_DPA1.3)
table(DPB1_fd$IMP_DPB1.3)
table(DRB1_fd$IMP_DRB1.3)
table(C_fd$IMP_C.3)
table(B_fd$IMP_B.3)
table(A_fd$IMP_A.3)
head(A_fd)


out <- merge(A_fd[,c('IID','IMP_A.1','IMP_A.2')],B_fd[,c('IID','IMP_B.1','IMP_B.2')],by = 'IID')
out <- merge(out,C_fd[,c('IID','IMP_C.1','IMP_C.2')],by = 'IID')
out <- merge(out,DRB1_fd[,c('IID','IMP_DRB1.1','IMP_DRB1.2')],by = 'IID')
out <- merge(out,DPA1_fd[,c('IID','IMP_DPA1.1','IMP_DPA1.2')],by = 'IID')
out <- merge(out,DPB1_fd[,c('IID','IMP_DPB1.1','IMP_DPB1.2')],by = 'IID')
out <- merge(out,DQA1_fd[,c('IID','IMP_DQA1.1','IMP_DQA1.2')],by = 'IID')
out <- merge(out,DQB1_fd[,c('IID','IMP_DQB1.1','IMP_DQB1.2')],by = 'IID')


write.csv(out,"HLAimputation.all.gene.4digit.Result.csv",row.names = F,quote = F)

head(out)
#####merge NGS
out <- read.csv("HLAimputation.all.gene.4digit.Result.csv")
ngs <- read.csv("../../../transplantation/HLAtyping/20200828/HLAtyping.alle.gene.4digit.csv")

head(ngs)
ncol(ngs)
out <- merge(out,ngs,by.x = "IID",by.y = "KID")
ncol(out)

head(out)
out <-out[,c(1,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]
write.csv(out,"MERGE.impResult.hlatyping.all.gene.4digit.csv",row.names = F,quote = F)








######python ÀÛ¾÷ ÈÄ..merge
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result3/ngs.vs.sm.compare/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result2/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/UsingHan/")
final_merge <- function(digit,A,B){
  a <- paste0("HLA.imputation.A_allele.",as.character(digit),"digit.without.3.allele.compare.",A,".and.",B,".txt")
  b <- paste0("HLA.imputation.B_allele.",as.character(digit),"digit.without.3.allele.compare.",A,".and.",B,".txt")
  drb <- paste0("HLA.imputation.DRB_allele.",as.character(digit),"digit.without.3.allele.compare.",A,".and.",B,".txt")
  a<-read.table(a,header = T)
  b<-read.table(b,header = T)
  drb<-read.table(drb,header = T)
  
  df <- merge(a,b,by = "ID")
  df <- merge(df,drb,by = "ID")
  colnames(df)[6:8] <-c("A.match","A.wrong","A.empty")
  colnames(df)[13:15] <-c("B.match","B.wrong","B.empty")
  colnames(df)[20:22] <-c("DRB1.match","DRB1.wrong","DRB1.empty")
  out = paste0("HLA.imputation.",as.character(digit),"digit.result.compare.",A,".and.",B,".csv")
  write.csv(df,out,row.names = F,quote = F)
  return(df)
}

df <- final_merge(2,"ngs","sm")
df <- final_merge(4,"ngs","sm")
ncol(df)
head(df)
