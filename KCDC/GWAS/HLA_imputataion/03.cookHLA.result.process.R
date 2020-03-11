setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result/")
df <- read.table("HLA_imputed_RAW.raw",header = T)
head(df)
colnames(df)
grep("*DPB1*",colnames(df))


library(stringr)

A <- df[,c(1,grep("*_A_*",colnames(df)))]
B <- df[,c(1,grep("*_B_*",colnames(df)))]
DPB <- df[,c(1,2,grep("*DPB1*",colnames(df)))]
head(DPB)


hla.subset <- function(a,n){
  if( (n ==2) || (n == 4)){
    i = paste(n,"(nvalue) is OK")
    print(i)
    temp <- a[,1:2]
    if(n == 2){
      for (i in 3:ncol(a)){
        if(nchar(colnames(a)[i]) == 13){
          b <- a[,c(1,2,i)]
          temp <- merge(temp,b)
        }
      }
    }else{
      for (i in 3:ncol(a)){
        if(nchar(colnames(a)[i]) > 13){
          b <- a[,c(1,2,i)]
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
DPB_td <- hla.subset(DPB,3)
DPB_td <- hla.subset(DPB,2)

head(DPB_td)



head(DPB_A2)


nchar(colnames(DPB)[-1]) == 13
DPB_A1 <- DPB[,c(1,(nchar(colnames(DPB)[-1]) == 13))]
head(DPB_A1)


  
a <- (nchar(colnames(DPB)[-1]) == 13)  
a <- "HLA_DPB1_0402_P"
str_split_fixed(a,"_",4)[,3]
str_split_fixed(a,"_",4)[,3]
nchar(a)








24000/4 + (16000 + 8900 + 43000)/5
(16000 + 8900 + 43000)/5

24000 - 19580 
16000 - 19580
