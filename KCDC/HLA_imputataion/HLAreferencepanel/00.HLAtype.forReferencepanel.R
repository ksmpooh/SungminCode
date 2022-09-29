library(tidyverse)
library(readxl)

setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/NGS(HLAtyping결과)/")
df <-read_excel("HLA.typing.Final.result_modify_20211216.xlsx")
head(df)
dim(df)


#FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1
#/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping/all/HLA.type.result.8genes.merged(2019.2020).csv

df <- read.csv("/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping/all/HLA.type.result.8genes.merged(2019.2020).csv")
head(df)
head(df[1:5,c("NGS_DPA1.12","NGS_DPA1.2")])

tcol <- grep("NGS_",colnames(df))

df[,tcol] <- lapply(df[tcol], function(x)
  str_split_fixed(x,"/",2)[,1]
  )
head(df)

df[ df == "." ] <- NA

df[(!is.na(df$NGS_DQA1.1)) & (is.na(df$NGS_DQA1.2)),]$NGS_DQA1.2 <- df[(!is.na(df$NGS_DQA1.1)) & (is.na(df$NGS_DQA1.2)),]$NGS_DQA1.1
df[(is.na(df$NGS_DQA1.1)) & (!is.na(df$NGS_DQA1.2)),]$NGS_DQA1.1 <- df[(is.na(df$NGS_DQA1.1)) & (!is.na(df$NGS_DQA1.2)),]$NGS_DQA1.2


df[(!is.na(df$NGS_DRB1.1)) & (is.na(df$NGS_DRB1.2)),]$NGS_DRB1.2 <- df[(!is.na(df$NGS_DRB1.1)) & (is.na(df$NGS_DRB1.2)),]$NGS_DRB1.1
df[(is.na(df$NGS_DRB1.1)) & (!is.na(df$NGS_DRB1.2)),]$NGS_DRB1.1 <- df[(is.na(df$NGS_DRB1.1)) & (!is.na(df$NGS_DRB1.2)),]$NGS_DRB1.2
#df[(!is.na(df$NGS_DRB3.1)) & (is.na(df$NGS_DRB3.2)),]$NGS_DRB3.2 <- df[(!is.na(df$NGS_DRB3.1)) & (is.na(df$NGS_DRB3.2)),]$NGS_DRB3.1
#df[(is.na(df$NGS_DRB3.1)) & (!is.na(df$NGS_DRB3.2)),]$NGS_DRB3.1 <- df[(is.na(df$NGS_DRB3.1)) & (!is.na(df$NGS_DRB3.2)),]$NGS_DRB3.2

df[(!is.na(df$NGS_DQA1.1)) & (is.na(df$NGS_DQA1.2)),]$NGS_DQA1.2 <- df[(!is.na(df$NGS_DQA1.1)) & (is.na(df$NGS_DQA1.2)),]$NGS_DQA1.1
df[(is.na(df$NGS_DQA1.1)) & (!is.na(df$NGS_DQA1.2)),]$NGS_DQA1.1 <- df[(is.na(df$NGS_DQA1.1)) & (!is.na(df$NGS_DQA1.2)),]$NGS_DQA1.2
df[(!is.na(df$NGS_DQB1.1)) & (is.na(df$NGS_DQB1.2)),]$NGS_DQB1.2 <- df[(!is.na(df$NGS_DQB1.1)) & (is.na(df$NGS_DQB1.2)),]$NGS_DQB1.1
df[(is.na(df$NGS_DQB1.1)) & (!is.na(df$NGS_DQB1.2)),]$NGS_DQB1.1 <- df[(is.na(df$NGS_DQB1.1)) & (!is.na(df$NGS_DQB1.2)),]$NGS_DQB1.2

df[(!is.na(df$NGS_DPA1.1)) & (is.na(df$NGS_DPA1.2)),]$NGS_DPA1.2 <- df[(!is.na(df$NGS_DPA1.1)) & (is.na(df$NGS_DPA1.2)),]$NGS_DPA1.1
df[(is.na(df$NGS_DPA1.1)) & (!is.na(df$NGS_DPA1.2)),]$NGS_DPA1.1 <- df[(is.na(df$NGS_DPA1.1)) & (!is.na(df$NGS_DPA1.2)),]$NGS_DPA1.2
df[(!is.na(df$NGS_DPB1.1)) & (is.na(df$NGS_DPB1.2)),]$NGS_DPB1.2 <- df[(!is.na(df$NGS_DPB1.1)) & (is.na(df$NGS_DPB1.2)),]$NGS_DPB1.1
df[(is.na(df$NGS_DPB1.1)) & (!is.na(df$NGS_DPB1.2)),]$NGS_DPB1.1 <- df[(is.na(df$NGS_DPB1.1)) & (!is.na(df$NGS_DPB1.2)),]$NGS_DPB1.2

df[(!is.na(df$NGS_A.1)) & (is.na(df$NGS_A.2)),]$NGS_A.2 <- df[(!is.na(df$NGS_A.1)) & (is.na(df$NGS_A.2)),]$NGS_A.1
df[(is.na(df$NGS_A.1)) & (!is.na(df$NGS_A.2)),]$NGS_A.1 <- df[(is.na(df$NGS_A.1)) & (!is.na(df$NGS_A.2)),]$NGS_A.2

df[(!is.na(df$NGS_B.1)) & (is.na(df$NGS_B.2)),]$NGS_B.2 <- df[(!is.na(df$NGS_B.1)) & (is.na(df$NGS_B.2)),]$NGS_B.1
df[(is.na(df$NGS_B.1)) & (!is.na(df$NGS_B.2)),]$NGS_B.1 <- df[(is.na(df$NGS_B.1)) & (!is.na(df$NGS_B.2)),]$NGS_B.2

df[(!is.na(df$NGS_C.1)) & (is.na(df$NGS_C.2)),]$NGS_C.2 <- df[(!is.na(df$NGS_C.1)) & (is.na(df$NGS_C.2)),]$NGS_C.1
df[(is.na(df$NGS_C.1)) & (!is.na(df$NGS_C.2)),]$NGS_C.1 <- df[(is.na(df$NGS_C.1)) & (!is.na(df$NGS_C.2)),]$NGS_C.2


#ngs.6d[(!is.na(ngs.6d$DRB1.1)) & (is.na(ngs.6d$DRB1.2)),]$DRB1.2 
#<- ngs.6d[(!is.na(ngs.6d$DRB1.1)) & (is.na(ngs.6d$DRB1.2)),]$DRB1.1


#FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1

df$FID <- df$KID
df$IID <- df$KID
df$pID <- 0
df$mID <- 0
df$SEX <- 0
df$PHENO <- 0
colnames(df)
# "NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2","NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DRB1.1","NGS_DRB1.2"
df <- df[,c("FID","IID","pID","mID","SEX","PHENO","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2","NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DRB1.1","NGS_DRB1.2")] 





#df <- read.csv("/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping/all/HLA.type.result.8genes.merged(2019.2020)_forMAKEreference.csv")
#write.table(df,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping/all/HLA.type.result.8genes.merged(2019.2020)_forMAKEreference.txt",col.names = T,row.names = F,quote = F)
#write.table(df,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLA.type.result.8genes.merged(2019.2020)_forMAKEreference.txt",col.names = T,row.names = F,quote = F,sep="\t")


#str_split_fixed(paste0(str_split_fixed(ngs$DPA1_1,":",3)[,1],str_split_fixed(ngs$DPA1_1,":",3)[,2]),"\\*",2)[,2]
head(df$NGS_A.1)

str_split_fixed(head(df$NGS_A.1),":",4)
paste0(str_split_fixed(head(df$NGS_A.1),":",3)[,1],":",str_split_fixed(head(df$NGS_A.1),":",3)[,2])
tcol <- grep("NGS_",colnames(df))
tcol
df[,tcol] <- lapply(df[tcol], function(x)
  paste0(str_split_fixed(x,":",3)[,1],":",str_split_fixed(x,":",3)[,2])
)
head(df)

write.table(df,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLA.type.result.8genes.merged.4digit_forMAKEreference.txt",col.names = T,row.names = F,quote = F,sep="\t")



tcol <- grep("NGS_",colnames(df))
tcol
df[,tcol] <- lapply(df[tcol], function(x)
  paste0(str_split_fixed(x,":",4)[,1],":",str_split_fixed(x,":",4)[,2],":",str_split_fixed(x,":",4)[,3])
)
head(df)

