### 특허 관련 hla imputation 결과 정리
setwd("~/Desktop/KCDC/HLAimputation/patent/split_result/")
####### Fianl allgene result processing using recodeA
gene = "A"
#gene = "B"
#gene = "DRB1"
theme = ""
#theme = "_pruning"
# mac

gene = "A"
#theme = ""
theme = "_pruning"
A <- read.table(paste0("JG.QCed.HLA_intersect_HLA.",gene,theme,"_raw.raw"),header = T)
header <- read.csv("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/Result/last.header.txt")
header <- header[,c(1,2,3,4,5,6,grep(paste0("HLA_",gene,"_*"),colnames(header)))]
colnames(A) <- colnames(header)
A[is.na(A)]<- -1
A <- A[,c(1:2,7:ncol(A))]



gene = "B"
#theme = ""
#theme = "_pruning"
B <- read.table(paste0("JG.QCed.HLA_intersect_HLA.",gene,theme,"_raw.raw"),header = T)
header <- read.csv("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/Result/last.header.txt")
header <- header[,c(1,2,3,4,5,6,grep(paste0("HLA_",gene,"_*"),colnames(header)))]
colnames(B) <- colnames(header)
B[is.na(B)]<- -1
B <- B[,c(1:2,7:ncol(B))]

gene = "DRB1"
#theme = ""
#theme = "_pruning"
DRB1 <- read.table(paste0("JG.QCed.HLA_intersect_HLA.",gene,theme,"_raw.raw"),header = T)
header <- read.csv("~/Desktop/KCDC/HLAimputation/IMPUTE4/Han.ref/Result/last.header.txt")
header <- header[,c(1,2,3,4,5,6,grep(paste0("HLA_",gene,"_*"),colnames(header)))]
colnames(DRB1) <- colnames(header)
DRB1[is.na(DRB1)]<- -1
DRB1 <- DRB1[,c(1:2,7:ncol(DRB1))]



library(stringr)

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
#### 2digit


###################### 2digit
###################### 2digit
grep("*N_P",colnames(A))
colnames(A)[6]<-"HLA_A_0122_P"
colnames(A)[19]<-"HLA_A_0253_P"  
####################
head(A)
A_td <- hla.subset(A,2)
A_td <- hla.find(A_td,"A")
head(A_td)

B_td <- hla.subset(B,2)
B_td <- hla.find(B_td,"B")
head(B_td)
DRB1_td <- hla.subset(DRB1,2)
DRB1_td <- hla.find(DRB1_td,"DRB1")
head(DRB1_td)




out <- merge(A_td[,c('IID','IMP_A.1','IMP_A.2')],B_td[,c('IID','IMP_B.1','IMP_B.2')],by = 'IID')
out <- merge(out,DRB1_td[,c('IID','IMP_DRB1.1','IMP_DRB1.2')],by = 'IID')

write.csv(out,"compare/HLA.gene.split.2d.csv",row.names = F,quote = F)

head(out)
#####merge NGS
ngs <- read.csv("../../HLAtyping/all/HLAtyping.alle.gene.2digit_2019.with.2020.csv")
head(ngs)
ncol(ngs)
out <- merge(out,ngs,by.x = "IID",by.y = "KID")
ncol(out)
head(out)

#out <-out[,c(1,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]
out <- out[,c("IID","IMP_A.1","IMP_A.2","NGS_A.1","NGS_A.2","IMP_B.1","IMP_B.2","NGS_B.1","NGS_B.2","IMP_DRB1.1","IMP_DRB1.2","NGS_DRB1.1","NGS_DRB1.2")]
#row.names(out) <- out$IID
#out[DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID,]
#out[DPA1_td[!is.na(DPA1_td$IMP_DPA1.3),]$IID,]
head(out)
nrow(out)
#write.csv(out,"compare/MERGE.impResult.hlatyping.splitimp.A.B.DRB1.gene.2digit.csv",row.names = F,quote = F)
write.csv(out,"compare/Pruning.MERGE.impResult.hlatyping.splitimp.A.B.DRB1.gene.2digit.csv",row.names = F,quote = F)






####################4 digit
A_fd <- hla.subset(A,4)
A_fd <- hla.find(A_fd,"A")
head(A_fd)
table(A_fd$IMP_A.3)
B_fd <- hla.subset(B,4)
B_fd <- hla.find(B_fd,"B")
head(B_fd)
DRB1_fd <- hla.subset(DRB1,4)
DRB1_fd <- hla.find(DRB1_fd,"DRB1")
head(DRB1_fd)

out <- merge(A_fd[,c('IID','IMP_A.1','IMP_A.2')],B_td[,c('IID','IMP_B.1','IMP_B.2')],by = 'IID')
out <- merge(out,DRB1_td[,c('IID','IMP_DRB1.1','IMP_DRB1.2')],by = 'IID')

write.csv(out,"compare/HLA.gene.split.4d.csv",row.names = F,quote = F)
write.csv(out,"compare/Pruning.HLA.gene.split.4d.csv",row.names = F,quote = F)


head(out)
#####merge NGS
ngs <- read.csv("../../HLAtyping/all/HLAtyping.alle.gene.4digit_2019.with.2020.csv")
head(ngs)
ncol(ngs)
out <- merge(out,ngs,by.x = "IID",by.y = "KID")
ncol(out)
head(out)

#out <-out[,c(1,18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]
out <- out[,c("IID","IMP_A.1","IMP_A.2","NGS_A.1","NGS_A.2","IMP_B.1","IMP_B.2","NGS_B.1","NGS_B.2","IMP_DRB1.1","IMP_DRB1.2","NGS_DRB1.1","NGS_DRB1.2")]
#row.names(out) <- out$IID
#out[DQA1_td[!is.na(DQA1_td$IMP_DQA1.3),]$IID,]
#out[DPA1_td[!is.na(DPA1_td$IMP_DPA1.3),]$IID,]
head(out)
nrow(out)
write.csv(out,"compare/MERGE.impResult.hlatyping.splitimp.A.B.DRB1.gene.4digit.csv",row.names = F,quote = F)
write.csv(out,"compare/Pruning.MERGE.impResult.hlatyping.splitimp.A.B.DRB1.gene.4digit.csv",row.names = F,quote = F)










