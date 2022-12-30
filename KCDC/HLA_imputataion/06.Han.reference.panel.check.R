setwd("c:/Users/user/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/Result/")
df <- read.table("han.impute4.hlaIMP_HLA_raw.raw",header = T)
header <- read.table("last.header.txt",sep = ",")
colnames(df)<-header
head(df)


###############
a <- table(is.na(df[,14]))
a <- as.matrix(table(is.na(df[,14])))
a <- as.data.frame(t(a))
a
rownames(a) <- NULL

colnames(a)<- c("NA","no.NA")
b <- data.frame(type = colnames(df)[14])
b
c <- cbind(b,a)
c
#############

out <- data.frame("type","NA","no.NA")
colnames(out) <- c("type","NA","no.NA")
out <- data.frame()


for (i in 7:ncol(df)) {
  a <- as.matrix(table(is.na(df[,i])))
  a <- as.data.frame(t(a))
  rownames(a) <- NULL
  b <- data.frame(type = colnames(df)[i])
  c <- cbind(b,a)
  out <- merge(out,c,all = TRUE)
}
out



#############################################################################

#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/HAN.ref/")
setwd("~/Desktop/KCDC/HLAimputation/HAN.ref/")
df <- read.table("HAN.hlatype.raw",header = T)
head(df)

library(stringr)


A <- df[,c(1,2,grep("HLA_A_*",colnames(df)))]
B <- df[,c(1,2,grep("HLA_B_*",colnames(df)))]
C <- df[,c(1,2,grep("HLA_C_*",colnames(df)))]
DRB1 <- df[,c(1,2,grep("*DRB1*",colnames(df)))]

DPB1 <- df[,c(1,2,grep("*DPB1*",colnames(df)))]
DPA1 <- df[,c(1,2,grep("*DPA1*",colnames(df)))]

DQB1 <- df[,c(1,2,grep("*DQB1*",colnames(df)))]
DQA1 <- df[,c(1,2,grep("*DQA1*",colnames(df)))]

head(A)

hla.subset <- function(df,n){
  if( (n ==2) || (n == 4)){
    i = paste(n,"(nvalue) is OK")
    print(i)
    temp <- df[,1:2]
    #    print(temp)
    if(n == 2){
      for (i in 3:ncol(df)){
        #print(as.integer(str_split_fixed(colnames(df)[i],"_",4)[3]))
        if(nchar(str_split_fixed(colnames(df)[i],"_",4)[3]) < 3){
          b <- df[,c(1,2,i)]
          temp <- merge(temp,b)
        }
      }
    }else{
      for (i in 3:ncol(df)){
        if(nchar(str_split_fixed(colnames(df)[i],"_",4)[3]) >= 3){
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
          df[i,paste("HLAtype_",concept,'.1',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          df[i,paste("HLAtype_",concept,'.2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 2  
        }else if(n == 1){
          df[i,paste("HLAtype_",concept,'.2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          df[i,paste("HLAtype_",concept,'.3',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 2  
        }
      }else if(df[i,j] == 1){
        if(n == 0){
          df[i,paste("HLAtype_",concept,'.1',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 1
        }else if(n == 1){
          df[i,paste("HLAtype_",concept,'.2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 1
        }else{
          df[i,paste("HLAtype_",concept,'.3',sep = "")]  <- str_split_fixed(colnames(df)[j],"_",4)[3]
        }
      }
    }
  }
  print(head(df))
  return(df)
}

DRB1_td <- hla.subset(DRB1,2)
DRB1_td <- hla.find(DRB1_td,"DRB1")
DQA1_td <- hla.subset(DQA1,2)
DQA1_td <- hla.find(DQA1_td,"DQA1")
DQB1_td <- hla.subset(DQB1,2)
DQB1_td <- hla.find(DQB1_td,"DQB1")

DPA1_td <- hla.subset(DPA1,2)
DPA1_td <- hla.find(DPA1_td,"DPA1")
DPB1_td <- hla.subset(DPB1,2)
DPB1_td <- hla.find(DPB1_td,"DPB1")

A_td <- hla.subset(A,2)
A_td <- hla.find(A_td,"A")
B_td <- hla.subset(B,2)
B_td <- hla.find(B_td,"B")
C_td <- hla.subset(C,2)
C_td <- hla.find(C_td,"C")

out <- merge(A_td[,c('IID','HLAtype_A.1','HLAtype_A.2')],B_td[,c('IID','HLAtype_B.1','HLAtype_B.2')],by = 'IID')
out <- merge(out,C_td[,c('IID','HLAtype_C.1','HLAtype_C.2')],by = 'IID')
out <- merge(out,DRB1_td[,c('IID','HLAtype_DRB1.1','HLAtype_DRB1.2')],by = 'IID')
out <- merge(out,DPA1_td[,c('IID','HLAtype_DPA1.1','HLAtype_DPA1.2')],by = 'IID')
out <- merge(out,DPB1_td[,c('IID','HLAtype_DPB1.1','HLAtype_DPB1.2')],by = 'IID')
out <- merge(out,DQA1_td[,c('IID','HLAtype_DQA1.1','HLAtype_DQA1.2')],by = 'IID')
out <- merge(out,DQB1_td[,c('IID','HLAtype_DQB1.1','HLAtype_DQB1.2')],by = 'IID')

####################################################################


DRB1_fd <- hla.subset(DRB1,4)
DRB1_fd <- hla.find(DRB1_fd,"DRB1")
DQA1_fd <- hla.subset(DQA1,4)
DQA1_fd <- hla.find(DQA1_fd,"DQA1")
DQB1_fd <- hla.subset(DQB1,4)
DQB1_fd <- hla.find(DQB1_fd,"DQB1")
DPA1_fd <- hla.subset(DPA1,4)
DPA1_fd <- hla.find(DPA1_fd,"DPA1")
DPB1_fd <- hla.subset(DPB1,4)
DPB1_fd <- hla.find(DPB1_fd,"DPB1")
A_fd <- hla.subset(A,4)
A_fd <- hla.find(A_fd,"A")
B_fd <- hla.subset(B,4)
B_fd <- hla.find(B_fd,"B")
C_fd <- hla.subset(C,4)
C_fd <- hla.find(C_fd,"C")

table(DQA1_fd$HLAtype_DQA1.3)
table(DQB1_fd$HLAtype_DQB1.3)
table(DPA1_fd$HLAtype_DPA1.3)
table(DPB1_fd$HLAtype_DPB1.3)
table(DRB1_fd$HLAtype_DRB1.3)
table(C_fd$HLAtype_C.3)
table(B_fd$HLAtype_B.3)
table(A_fd$HLAtype_A.3)

td <- out
#head(out)

out <- merge(A_fd[,c('IID','HLAtype_A.1','HLAtype_A.2')],B_fd[,c('IID','HLAtype_B.1','HLAtype_B.2')],by = 'IID')
out <- merge(out,C_fd[,c('IID','HLAtype_C.1','HLAtype_C.2')],by = 'IID')
out <- merge(out,DRB1_fd[,c('IID','HLAtype_DRB1.1','HLAtype_DRB1.2')],by = 'IID')
out <- merge(out,DPA1_fd[,c('IID','HLAtype_DPA1.1','HLAtype_DPA1.2')],by = 'IID')
out <- merge(out,DPB1_fd[,c('IID','HLAtype_DPB1.1','HLAtype_DPB1.2')],by = 'IID')
out <- merge(out,DQA1_fd[,c('IID','HLAtype_DQA1.1','HLAtype_DQA1.2')],by = 'IID')
out <- merge(out,DQB1_fd[,c('IID','HLAtype_DQB1.1','HLAtype_DQB1.2')],by = 'IID')
fd <- out


head(td)

table(td$HLAtype_A.1)
a <- as.data.frame(table(td$HLAtype_A.1))
a
b <- as.data.frame(table(td$HLAtype_A.2))
b

df <- merge(a,b,all = T)
head(df)
table(td$HLAtype_A.1) + table(td$HLAtype_A.2)

head(fd)

write.csv(td,"Han.HLAtyping.Result.2digit.csv",col.names = T,row.names = F,quote = F)
write.csv(fd,"Han.HLAtyping.Result.4digit.csv",col.names = T,row.names = F,quote = F)



#####
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/HAN.ref/")
td <- read.csv("Han.HLAtyping.Result.2digit.csv")
fd <- read.csv("Han.HLAtyping.Result.4digit.csv")
head(td)
table(td$HLAtype_A.1)
d