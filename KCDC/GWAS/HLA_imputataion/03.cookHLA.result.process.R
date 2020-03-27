setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result3/")
df <- read.table("JG.HLA.imputation_RAW.raw",header = T)
head(df)
colnames(df)
grep("*DRB1*",colnames(df))


library(stringr)

A <- df[,c(1,2,grep("*_A_*",colnames(df)))]
B <- df[,c(1,2,grep("*_B_*",colnames(df)))]
DRB <- df[,c(1,2,grep("*DRB1*",colnames(df)))]
head(DRB)

hla.subset <- function(df,n){
  if( (n ==2) || (n == 4)){
    i = paste(n,"(nvalue) is OK")
    print(i)
    temp <- df[,1:2]
#    print(temp)
    if(n == 2){
      for (i in 3:ncol(df)){
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
####
a <- DRB_td[!is.na(DRB_td$DRB3),]
write.table(a,"HLAimputation.3alleles.in.DRB.2digit.txt",col.names = T,row.names = F,quote = F)
#####

DRB_td <- DRB_td[is.na(DRB_td$DRB3),]

DRB_td <- DRB_td[,c(1,2,ncol(DRB_td)-2,ncol(DRB_td)-1,ncol(DRB_td))]
#DRB_td <- DRB_td[,c(1,2,ncol(DRB_td)-2,ncol(DRB_td)-1)]
head(DRB_td)
table(DRB_td$DRB3)


A_td <- hla.subset(A,2)
A_td <- hla.find(A_td,"A")

A_td <- A_td[,c(1,2,ncol(A_td)-1,ncol(A_td))]
head(A_td)


B_td <- hla.subset(B,2)
B_td <- hla.find(B_td,"B")
#######
b <- B_td[!is.na(B_td$B3),]
table(b$B3)
write.table(b,"HLAimputation.3alleles.in.B.2digit.txt",col.names = T,row.names = F,quote = F)

####
B_td <- B_td[is.na(B_td$B3),]
B_td <- B_td[,c(1,2,ncol(B_td)-2,ncol(B_td)-1,ncol(B_td))]
#B_td <- B_td[,c(1,2,ncol(B_td)-2,ncol(B_td)-1)]

table(B_td$B3)
head(B_td)
####################################################################################################

df <- merge(A_td,B_td)
df <- merge(df,DRB_td)
head(df)

df <-df[is.na(df$B3),]
df <-df[is.na(df$DRB3),]
df <- df[,c(2,3,4,5,6,8,9)]
colnames(df)[1] <- "ID"

df[,2:ncol(df)] <- sapply(df[,2:ncol(df)],as.integer)
str(df)

ref <- read.table("../../transplantation/HLAtyping/HLA_JG_2DGT_imputed.txt",header = T)

a <- merge(df[,c(1,2,3)],ref[,c(1,2,3)],by = "ID",all.x = T)
b <- merge(df[,c(1,4,5)],ref[,c(1,4,5)],by = "ID",all.x = T)
drb <- merge(df[,c(1,6,7)],ref[,c(1,6,7)],by = "ID",all.x = T)

colnames(a) <- c("ID","sm.A1","sm.A2","yj.A1","yj.A2")
colnames(b) <- c("ID","sm.B1","sm.B2","yj.B1","yj.B2")
colnames(drb) <- c("ID","sm.DRB1","sm.DRB2","yj.DRB1","yj.DRB2")
head(a)
dim(a)
dim(b)
dim(drb)

write.table(a,"HLA_imputation_A_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')
write.table(b,"HLA_imputation_B_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')
write.table(drb,"HLA_imputation_DRB_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')

#write.table(df[,c(2,3,4,5,6,8,9)],"HLA_imptation_2d_without_3_alleles.txt",col.names = T,row.names = F,quote = F)



#############################################################################################


DRB_fd <- hla.subset(DRB,4)
head(DRB_fd)
DRB_fd <- hla.find(DRB_fd,"DRB")


DRB_fd <- DRB_fd[,c(1,2,ncol(DRB_fd)-2,ncol(DRB_fd)-1,ncol(DRB_fd))]
head(DRB_fd)
table(DRB_fd$DRB3)

a <- DRB_fd[!is.na(DRB_fd$DRB3),]
write.table(a,"HLAimputation.3alleles.in.DRB.4digit.txt",col.names = T,row.names = F,quote = F)

#DRB_fd <- DRB_fd[is.na(DRB_fd$DRB3),]

A_fd <- hla.subset(A,4)
A_fd <- hla.find(A_fd,"A")
A_fd <- A_fd[,c(1,2,ncol(A_fd)-1,ncol(A_fd))]
head(A_fd)


B_fd <- hla.subset(B,4)
B_fd <- hla.find(B_fd,"B")
head(B_fd)
B_fd <- B_fd[,c(1,2,ncol(B_fd)-1,ncol(B_fd))]

df <- merge(A_fd,B_fd)
df <- merge(df,DRB_fd)


df <-df[is.na(df$DRB3),]
head(df)
df <- subset(df,select = -c(DRB3,FID))

colnames(df)[1] <- "ID"

df[,2:ncol(df)] <- sapply(df[,2:ncol(df)],as.integer)
str(df)

#ref <- read.table("../../transplantation/HLAtyping/HLA_JG_4DGT_imputed.txt",header = T)
ref <- read.csv("../../transplantation/HLAtyping/HLA_JG_4DGT_imputed.csv",header = T)
head(ref)

a <- merge(df[,c(1,2,3)],ref[,c(1,2,3)],by = "ID",all.x = T)
b <- merge(df[,c(1,4,5)],ref[,c(1,4,5)],by = "ID",all.x = T)
drb <- merge(df[,c(1,6,7)],ref[,c(1,6,7)],by = "ID",all.x = T)

colnames(a) <- c("ID","sm.A1","sm.A2","yj.A1","yj.A2")
colnames(b) <- c("ID","sm.B1","sm.B2","yj.B1","yj.B2")
colnames(drb) <- c("ID","sm.DRB1","sm.DRB2","yj.DRB1","yj.DRB2")
head(a)
dim(a)
dim(b)
dim(drb)

write.table(a,"HLA_imputation_A_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')
write.table(b,"HLA_imputation_B_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')
write.table(drb,"HLA_imputation_DRB_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = '\t')



##################python data processing 후 merge 작업
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result3/")
final_merge <- function(digit){
  
  a <- paste0("HLA.imputation.A_allele.",as.character(digit),"digit.without.3.allele.compare.sm.and.yj.txt")
  b <- paste0("HLA.imputation.B_allele.",as.character(digit),"digit.without.3.allele.compare.sm.and.yj.txt")
  drb <- paste0("HLA.imputation.DRB_allele.",as.character(digit),"digit.without.3.allele.compare.sm.and.yj.txt")
  
  a<-read.table(a,header = T)
  b<-read.table(b,header = T)
  drb<-read.table(drb,header = T)
  
  df <- merge(a,b,by = "ID")
  df <- merge(df,drb,by = "ID")
  out = paste0("HLA.imputation.",as.character(digit),"digit.result.compare.sm.with.yj.csv")
  write.csv(df,out,row.names = F,quote = F)
  return(df)
}
df <- final_merge(2)
df <- final_merge(4)
head(df)













#####################처음에 했던 merge
a <- read.table("HLA.imputation.A_allele.2digit.without.3.allele.compare.sm.and.yj.txt",header = T)
b <- read.table("HLA.imputation.B_allele.2digit.without.3.allele.compare.sm.and.yj.txt",header = T)
drb <- read.table("HLA.imputation.DRB_allele.2digit.without.3.allele.compare.sm.and.yj.txt",header = T)
head(a)
head(b)
head(drb)


df <- merge(a,b,by = "ID")
df <- merge(df,drb,by = "ID")
head(df)
write.csv(df,"HLA.imputation.2digit.result.compare.sm.with.yj.csv",col.names = T,row.names = F,quote = F)

nrow(df)
sum(df$match.x + df$match + df$match.y)
sum(df$wrong.x + df$wrong.y + df$wrong)
sum(df$empty, df$empty.x, df$empty.y)

sum(df$match.x + df$match + df$match.y)/(sum(df$match.x + df$match + df$match.y) + sum(df$wrong.x + df$wrong.y + df$wrong))

  