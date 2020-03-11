setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result/")
df <- read.table("HLA_imputed_RAW.raw",header = T)
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
DRB_td <- hla.find(DRB_td,"DRB")
DRB_td <- DRB_td[,c(1,2,ncol(DRB_td)-2,ncol(DRB_td)-1,ncol(DRB_td))]
head(DRB_td)
table(DRB_td$DRB3)

A_td <- hla.subset(A,2)
A_td <- hla.find(A_td,"A")
A_td <- A_td[,c(1,2,ncol(A_td)-1,ncol(A_td))]
head(A_td)


B_td <- hla.subset(B,2)
B_td <- hla.find(B_td,"B")
B_td <- B_td[,c(1,2,ncol(B_td)-2,ncol(B_td)-1,ncol(B_td))]

table(B_td$B3)
head(B_td)
####################################################################################################

df <- merge(A_td,B_td)
df <- merge(df,DRB_td)
head(df)

write.table(df,"HLA_imptation_2d.txt",col.names = T,row.names = F,quote = F)



