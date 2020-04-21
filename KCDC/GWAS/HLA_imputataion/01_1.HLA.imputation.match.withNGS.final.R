#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/UsingHan/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result2/")

library(stringr)
ori <- read.csv("../../transplantation/HLAtyping/HLA_NGS_typing_255samples_results_202002.csv")

head(ori)               
colnames(ori)
head(ori)

ngs <-ori
head(ngs)
head(ngs)

paste0(str_split_fixed(ngs$A.Allele1,":",3)[,1],str_split_fixed(ngs$A.Allele1,":",3)[,2])
str_split_fixed(paste0(str_split_fixed(ngs$A.Allele1,":",3)[,1],str_split_fixed(ngs$A.Allele1,":",3)[,2]),"\\*",2)[,2]


##### 2digit..아래도있음
ngs$A1 <- str_split_fixed(str_split_fixed(ngs$A.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$A2 <- str_split_fixed(str_split_fixed(ngs$A.Allele2,":",2)[,1],"\\*",2)[,2]
ngs$B1 <- str_split_fixed(str_split_fixed(ngs$B.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$B2 <- str_split_fixed(str_split_fixed(ngs$B.Allele2,":",2)[,1],"\\*",2)[,2]
ngs$DRB1 <- str_split_fixed(str_split_fixed(ngs$DRB1.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$DRB2 <- str_split_fixed(str_split_fixed(ngs$DRB1.Allele2,":",2)[,1],"\\*",2)[,2]


ngs.2d <- subset(ngs,select =c("KID","A1","A2","B1","B2","DRB1","DRB2"))
features <- colnames(ngs.2d)[2:ncol(ngs.2d)]
features
for (i in features){
  #ngs.2d[,i] <- strtoi(ngs.2d[,i])
  ngs.2d[,i] <- as.integer(ngs.2d[,i])
}

ngs.2d[(!is.na(ngs.2d$DRB1)) & (is.na(ngs.2d$DRB2)),]$DRB2 <- ngs.2d[(!is.na(ngs.2d$DRB1)) & (is.na(ngs.2d$DRB2)),]$DRB1
ngs.2d[(is.na(ngs.2d$DRB1)) & (!is.na(ngs.2d$DRB2)),]$DRB1 <- ngs.2d[(is.na(ngs.2d$DRB1)) & (!is.na(ngs.2d$DRB2)),]$DRB2
ngs.2d[(!is.na(ngs.2d$A1)) & (is.na(ngs.2d$A2)),]$A2 <- ngs.2d[(!is.na(ngs.2d$A1)) & (is.na(ngs.2d$A2)),]$A1
ngs.2d[(is.na(ngs.2d$A1)) & (!is.na(ngs.2d$A2)),]$A1 <- ngs.2d[(is.na(ngs.2d$A1)) & (!is.na(ngs.2d$A2)),]$A2
ngs.2d[(!is.na(ngs.2d$B1)) & (is.na(ngs.2d$B2)),]$B2 <- ngs.2d[(!is.na(ngs.2d$B1)) & (is.na(ngs.2d$B2)),]$B1
ngs.2d[(is.na(ngs.2d$B1)) & (!is.na(ngs.2d$B2)),]$B1 <- ngs.2d[(is.na(ngs.2d$B1)) & (!is.na(ngs.2d$B2)),]$B2

#a[(!is.na(a$DRB1)) & (is.na(a$DRB2)),]$DRB2 <- a[(!is.na(a$DRB1)) & (is.na(a$DRB2)),]$DRB1
colnames(ngs.2d) <-c("ID","ngs.A1","ngs.A2","ngs.B1","ngs.B2","ngs.DRB1","ngs.DRB2")


A.2d <- read.table("HLA_imputation_A_allele_2d_without.3.allele.txt",header = T)
A.2d <- merge(ngs.2d[1:255,c(1,2,3)],A.2d[,c(1,2,3)],by = "ID",all.x = T)
head(A.2d)

write.table(A.2d,"ngs.vs.sm.compare/HLA_imputation_A_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")


B.2d <- read.table("HLA_imputation_B_allele_2d_without.3.allele.txt",header = T)
B.2d <- merge(ngs.2d[1:255,c(1,4,5)],B.2d[,c(1,2,3)],by = "ID",all.x = T)
head(B.2d)

write.table(B.2d,"ngs.vs.sm.compare/HLA_imputation_B_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")


DRB.2d <- read.table("HLA_imputation_DRB_allele_2d_without.3.allele.txt",header = T)
DRB.2d <- merge(ngs.2d[1:255,c(1,6,7)],DRB.2d[,c(1,2,3)],by = "ID",all.x = T)
head(DRB.2d)

write.table(DRB.2d,"ngs.vs.sm.compare/HLA_imputation_DRB_allele_2d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")



################4digit



##### 4digit
ngs$A1 <- str_split_fixed(paste0(str_split_fixed(ngs$A.Allele1,":",3)[,1],str_split_fixed(ngs$A.Allele1,":",3)[,2]),"\\*",2)[,2]
ngs$A2 <- str_split_fixed(paste0(str_split_fixed(ngs$A.Allele2,":",3)[,1],str_split_fixed(ngs$A.Allele2,":",3)[,2]),"\\*",2)[,2]
ngs$B1 <- str_split_fixed(paste0(str_split_fixed(ngs$B.Allele1,":",3)[,1],str_split_fixed(ngs$B.Allele1,":",3)[,2]),"\\*",2)[,2]
ngs$B2 <- str_split_fixed(paste0(str_split_fixed(ngs$B.Allele2,":",3)[,1],str_split_fixed(ngs$B.Allele2,":",3)[,2]),"\\*",2)[,2]
ngs$DRB1 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB1.Allele1,":",3)[,1],str_split_fixed(ngs$DRB1.Allele1,":",3)[,2]),"\\*",2)[,2]
ngs$DRB2 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB1.Allele2,":",3)[,1],str_split_fixed(ngs$DRB1.Allele2,":",3)[,2]),"\\*",2)[,2]


ngs.2d <- subset(ngs,select =c("KID","A1","A2","B1","B2","DRB1","DRB2"))
features <- colnames(ngs.2d)[2:ncol(ngs.2d)]
features
for (i in features){
  #ngs.2d[,i] <- strtoi(ngs.2d[,i])
  ngs.2d[,i] <- as.integer(ngs.2d[,i])
}
head(ngs.2d)

ngs.2d[(!is.na(ngs.2d$DRB1)) & (is.na(ngs.2d$DRB2)),]$DRB2 <- ngs.2d[(!is.na(ngs.2d$DRB1)) & (is.na(ngs.2d$DRB2)),]$DRB1
ngs.2d[(is.na(ngs.2d$DRB1)) & (!is.na(ngs.2d$DRB2)),]$DRB1 <- ngs.2d[(is.na(ngs.2d$DRB1)) & (!is.na(ngs.2d$DRB2)),]$DRB2
ngs.2d[(!is.na(ngs.2d$A1)) & (is.na(ngs.2d$A2)),]$A2 <- ngs.2d[(!is.na(ngs.2d$A1)) & (is.na(ngs.2d$A2)),]$A1
ngs.2d[(is.na(ngs.2d$A1)) & (!is.na(ngs.2d$A2)),]$A1 <- ngs.2d[(is.na(ngs.2d$A1)) & (!is.na(ngs.2d$A2)),]$A2

ngs.2d[(!is.na(ngs.2d$B1)) & (is.na(ngs.2d$B2)),]$B2 <- ngs.2d[(!is.na(ngs.2d$B1)) & (is.na(ngs.2d$B2)),]$B1
ngs.2d[(is.na(ngs.2d$B1)) & (!is.na(ngs.2d$B2)),]$B1 <- ngs.2d[(is.na(ngs.2d$B1)) & (!is.na(ngs.2d$B2)),]$B2

colnames(ngs.2d) <-c("ID","ngs.A1","ngs.A2","ngs.B1","ngs.B2","ngs.DRB1","ngs.DRB2")
ngs.2d[1:10,]


A.2d <- read.table("HLA_imputation_A_allele_4d_without.3.allele.txt",header = T)
A.2d <- merge(ngs.2d[1:255,c(1,2,3)],A.2d[,c(1,2,3)],by = "ID",all.x = T)
head(A.2d)

write.table(A.2d,"ngs.vs.sm.compare/HLA_imputation_A_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")


B.2d <- read.table("HLA_imputation_B_allele_4d_without.3.allele.txt",header = T)
B.2d <- merge(ngs.2d[1:255,c(1,4,5)],B.2d[,c(1,2,3)],by = "ID",all.x = T)
head(B.2d)

write.table(B.2d,"ngs.vs.sm.compare/HLA_imputation_B_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")


DRB.2d <- read.table("HLA_imputation_DRB_allele_4d_without.3.allele.txt",header = T)
DRB.2d <- merge(ngs.2d[1:255,c(1,6,7)],DRB.2d[,c(1,2,3)],by = "ID",all.x = T)
head(DRB.2d)

write.table(DRB.2d,"ngs.vs.sm.compare/HLA_imputation_DRB_allele_4d_without.3.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")








######python 작업 후..merge
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result3/ngs.vs.sm.compare/")
final_merge <- function(digit,A,B){
  a <- paste0("ngs.vs.sm.compare/HLA.imputation.A_allele.",as.character(digit),"digit.without.3.allele.compare.",A,".and.",B,".txt")
  b <- paste0("ngs.vs.sm.compare/HLA.imputation.B_allele.",as.character(digit),"digit.without.3.allele.compare.",A,".and.",B,".txt")
  drb <- paste0("ngs.vs.sm.compare/HLA.imputation.DRB_allele.",as.character(digit),"digit.without.3.allele.compare.",A,".and.",B,".txt")
  a<-read.table(a,header = T)
  b<-read.table(b,header = T)
  drb<-read.table(drb,header = T)
  
  df <- merge(a,b,by = "ID")
  df <- merge(df,drb,by = "ID")
  colnames(df)[6:8] <-c("A.match","A.wrong","A.empty")
  colnames(df)[13:15] <-c("B.match","B.wrong","B.empty")
  colnames(df)[20:22] <-c("DRB.match","DRB.wrong","DRB.empty")
  out = paste0("ngs.vs.sm.compare/HLA.imputation.",as.character(digit),"digit.result.compare.",A,".and.",B,".csv")
  write.csv(df,out,row.names = F,quote = F)
  return(df)
}

df <- final_merge(2,"ngs","sm")
df <- final_merge(4,"ngs","sm")
ncol(df)
head(df)
