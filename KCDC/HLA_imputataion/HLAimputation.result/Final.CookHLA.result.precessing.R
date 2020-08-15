### FINAL cookHLA result processing and compare with NGS

#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result3/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/UsingHan/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20200731/Han/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/20200731/Pan/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/01.pan/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/255sample/02.han/")
#df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
head(df)
colnames(df)
#grep("*DRB1*",colnames(df))


library(stringr)

A <- df[,c(1,2,grep("HLA_A_*",colnames(df)))]
B <- df[,c(1,2,grep("HLA_B_*",colnames(df)))]

DRB <- df[,c(1,2,grep("*DRB1*",colnames(df)))]
head(DRB)
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
          df[i,paste(concept,'1',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          df[i,paste(concept,'2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 2  
        }else if(n == 1){
          df[i,paste(concept,'2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          df[i,paste(concept,'3',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 2  
        }
      }else if(df[i,j] == 1){
        if(n == 0){
          df[i,paste(concept,'1',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 1
        }else if(n == 1){
          df[i,paste(concept,'2',sep = "")] <- str_split_fixed(colnames(df)[j],"_",4)[3]
          n = n + 1
        }else{
          df[i,paste(concept,'3',sep = "")]  <- str_split_fixed(colnames(df)[j],"_",4)[3]
        }
      }
    }
  }
  print(head(df))
  return(df)
}

#hla.subset(df,digit)
#hla.find(df_digit,"concept")
DRB_td <- hla.subset(DRB,2)
head(DRB_td)
DRB_td <- hla.find(DRB_td,"DRB")
table(DRB_td$DRB3)
####DRB3이 있을경우
a <- DRB_td[!is.na(DRB_td$DRB3),]
write.table(a,"HLAimputation.3alleles.in.DRB.2digit.txt",col.names = T,row.names = F,quote = F)
DRB_td <- DRB_td[is.na(DRB_td$DRB3),]
#####



#DRB_td <- DRB_td[,c(1,2,ncol(DRB_td)-2,ncol(DRB_td)-1,ncol(DRB_td))]
#DRB_td <- DRB_td[,c(1,2,ncol(DRB_td)-2,ncol(DRB_td)-1)]
DRB_td <- subset(DRB_td,select = c("IID","DRB1","DRB2"))
head(DRB_td)
table(DRB_td$DRB3)


####20200417 usnig Han processing
grep("*N_P",colnames(A))
colnames(A)[6]<-"HLA_A_0122_P"
colnames(A)[19]<-"HLA_A_0253_P"

#colnames(A)[6]<-"HLA_A_0122_P"
colnames(A)[19]
###20200417

A_td <- hla.subset(A,2)
A_td <- hla.find(A_td,"A")
table(A_td$A3)
###A3가 있을 경우

a <- A_td[!is.na(A_td$A3),]
write.table(a,"HLAimputation.3alleles.in.A.2digit.txt",col.names = T,row.names = F,quote = F)
A_td <- A_td[is.na(A_td$A3),]


###
#A_td <- A_td[,c(1,2,ncol(A_td)-1,ncol(A_td))]
A_td <- A_td[,c("IID","A1","A2")]
head(A_td)
###
B_td <- hla.subset(B,2)
B_td <- hla.find(B_td,"B")
table(B_td$B3)
head(B_td)
#####B3가 있을경우
b <- B_td[!is.na(B_td$B3),]
table(b$B3)
write.table(b,"HLAimputation.3alleles.in.B.2digit.txt",col.names = T,row.names = F,quote = F)
B_td <- B_td[is.na(B_td$B3),]
####

#B_td <- B_td[,c(1,2,ncol(B_td)-2,ncol(B_td)-1,ncol(B_td))]
#B_td <- B_td[,c(1,2,ncol(B_td)-2,ncol(B_td)-1)]
B_td <- subset(B_td,select = c("IID","B1","B2"))
table(B_td$B3)
head(B_td)
####################################################################################################

######allele 3이 있을경우
#df <-df[is.na(df$B3),]
#df <-df[is.na(df$DRB3),]
#df <- df[,c(2,3,4,5,6,8,9)]
#df <- subset(df,select = c("IID","A1","A2","B1","B2","DRB1","DRB2"))

######2 allele 일경우


colnames(A_td) <- c("ID","sm.A1","sm.A2")
colnames(B_td) <- c("ID","sm.B1","sm.B2")
colnames(DRB_td) <- c("ID","sm.DRB1","sm.DRB2")
#colnames(df)[1] <- "ID"

head(A_td)

A_td[,2:ncol(A_td)] <- sapply(A_td[,2:ncol(A_td)],as.integer)
B_td[,2:ncol(B_td)] <- sapply(B_td[,2:ncol(B_td)],as.integer)
DRB_td[,2:ncol(DRB_td)] <- sapply(DRB_td[,2:ncol(DRB_td)],as.integer)


#write.table(A_td,"HLA_imputation_A_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')
#write.table(B_td,"HLA_imputation_B_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')
#write.table(DRB_td,"HLA_imputation_DRB_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')


ref <- read.csv("../../NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype_processing.2digit.csv",header = T)
head(ref)

A.2d <- merge(ref[,c("ID","ngs.A1","ngs.A2")],A_td[,c("ID","sm.A1","sm.A2")],by = "ID",all.x = T)
B.2d <- merge(ref[,c("ID","ngs.B1","ngs.B2")],B_td[,c("ID","sm.B1","sm.B2")],by = "ID",all.x = T)
DRB.2d <- merge(ref[,c("ID","ngs.DRB1","ngs.DRB2")],DRB_td[,c("ID","sm.DRB1","sm.DRB2")],by = "ID",all.x = T)
head(B.2d)
head(B_td)
write.table(A.2d,"HLA_imputation_A_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(B.2d,"HLA_imputation_B_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(DRB.2d,"HLA_imputation_DRB_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")


#############################################################################################


DRB_fd <- hla.subset(DRB,4)
head(DRB_fd)
DRB_fd <- hla.find(DRB_fd,"DRB")
table(DRB_fd$DRB3)

##### 3 allele
a <- DRB_fd[!is.na(DRB_fd$DRB3),]
write.table(a,"HLAimputation.3alleles.in.DRB.4digit.txt",col.names = T,row.names = F,quote = F)
DRB_fd <- DRB_fd[is.na(DRB_fd$DRB3),]
#####

DRB_fd <- DRB_fd[,c("FID","DRB1","DRB2")]
#DRB_fd <- DRB_fd[is.na(DRB_fd$DRB3),]

A_fd <- hla.subset(A,4)
A_fd <- hla.find(A_fd,"A")

head(A_fd$A3)
table(A_fd$A3)
##### 3 allele
A_fd <- A_fd[is.na(A_fd$A3),]
#######

A_fd <- A_fd[,c("FID","A1","A2")]




B_fd <- hla.subset(B,4)
B_fd <- hla.find(B_fd,"B")
head(B_fd)

##### 3 allele
B_fd <- B_fd[is.na(B_fd$A3),]
#####
B_fd <- subset(B_fd,select = c("FID","B1","B2"))



colnames(A_fd) <- c("ID","sm.A1","sm.A2")
colnames(B_fd) <- c("ID","sm.B1","sm.B2")
colnames(DRB_fd) <- c("ID","sm.DRB1","sm.DRB2")

head(A_fd)

A_fd[,2:ncol(A_fd)] <- sapply(A_fd[,2:ncol(A_fd)],as.integer)
B_fd[,2:ncol(B_fd)] <- sapply(B_fd[,2:ncol(B_fd)],as.integer)
DRB_fd[,2:ncol(DRB_fd)] <- sapply(DRB_fd[,2:ncol(DRB_fd)],as.integer)

#####################################################3

ref <- read.csv("../../NGS/HLA_NGS_typing_255samples_results_202002_modify_alltype_processing.4digit.csv",header = T)
head(ref)

A.4d <- merge(ref[,c("ID","ngs.A1","ngs.A2")],A_fd[,c("ID","sm.A1","sm.A2")],by = "ID",all.x = T)
B.4d <- merge(ref[,c("ID","ngs.B1","ngs.B2")],B_fd[,c("ID","sm.B1","sm.B2")],by = "ID",all.x = T)
DRB.4d <- merge(ref[,c("ID","ngs.DRB1","ngs.DRB2")],DRB_fd[,c("ID","sm.DRB1","sm.DRB2")],by = "ID",all.x = T)
head(B.4d)
head(B_fd)
write.table(A.4d,"HLA_imputation_A_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(B.4d,"HLA_imputation_B_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(DRB.4d,"HLA_imputation_DRB_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")






##################python data processing 후 merge 작업













######python 작업 후..merge
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
