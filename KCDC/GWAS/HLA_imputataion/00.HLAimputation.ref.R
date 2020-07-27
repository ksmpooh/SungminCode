
library(stringr)
ori <- read.csv("HLA_NGS_typing_result.csv")

##### 2digit..
ngs$A1 <- str_split_fixed(str_split_fixed(ngs$A.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$A2 <- str_split_fixed(str_split_fixed(ngs$A.Allele2,":",2)[,1],"\\*",2)[,2]
ngs$B1 <- str_split_fixed(str_split_fixed(ngs$B.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$B2 <- str_split_fixed(str_split_fixed(ngs$B.Allele2,":",2)[,1],"\\*",2)[,2]
ngs$DRB1 <- str_split_fixed(str_split_fixed(ngs$DRB1.Allele1,":",2)[,1],"\\*",2)[,2]
ngs$DRB2 <- str_split_fixed(str_split_fixed(ngs$DRB1.Allele2,":",2)[,1],"\\*",2)[,2]

#### 원하는 column 선택
ngs.2d <- subset(ngs,select =c("KID","A1","A2","B1","B2","DRB1","DRB2"))
features <- colnames(ngs.2d)[2:ncol(ngs.2d)]

#### string으로 되어 있는 value들을 integer로 변환
for (i in features){
  ngs.2d[,i] <- as.integer(ngs.2d[,i])
}

## homozygote allele 처리
ngs.2d[(!is.na(ngs.2d$DRB1)) & (is.na(ngs.2d$DRB2)),]$DRB2 <- ngs.2d[(!is.na(ngs.2d$DRB1)) & (is.na(ngs.2d$DRB2)),]$DRB1
ngs.2d[(is.na(ngs.2d$DRB1)) & (!is.na(ngs.2d$DRB2)),]$DRB1 <- ngs.2d[(is.na(ngs.2d$DRB1)) & (!is.na(ngs.2d$DRB2)),]$DRB2
ngs.2d[(!is.na(ngs.2d$A1)) & (is.na(ngs.2d$A2)),]$A2 <- ngs.2d[(!is.na(ngs.2d$A1)) & (is.na(ngs.2d$A2)),]$A1
ngs.2d[(is.na(ngs.2d$A1)) & (!is.na(ngs.2d$A2)),]$A1 <- ngs.2d[(is.na(ngs.2d$A1)) & (!is.na(ngs.2d$A2)),]$A2
ngs.2d[(!is.na(ngs.2d$B1)) & (is.na(ngs.2d$B2)),]$B2 <- ngs.2d[(!is.na(ngs.2d$B1)) & (is.na(ngs.2d$B2)),]$B1
ngs.2d[(is.na(ngs.2d$B1)) & (!is.na(ngs.2d$B2)),]$B1 <- ngs.2d[(is.na(ngs.2d$B1)) & (!is.na(ngs.2d$B2)),]$B2

## 추후 HLA imputation 결과값과 비교하기 위해 columne name 변경 
colnames(ngs.2d) <-c("ID","ngs.A1","ngs.A2","ngs.B1","ngs.B2","ngs.DRB1","ngs.DRB2")

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

######python 작업 후..merge
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result3/ngs.vs.sm.compare/")
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/Result2/")
#setwd("c:/Users/user/Desktop/KCDC/HLAimputation/UsingHan/")
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




######examplte

x = c("A:00*01*02*03","A:11*01*02*03")
str_split_fixed(x,":",2)
str_split_fixed(x,"\\*",2)
str_split_fixed(x,"\\*",4)
str_split_fixed(x,":",3)[,1]
str_split_fixed(x,":",3)[,2]
str_split_fixed(str_split_fixed(x,":",2)[,2],"\\*",2)[,1]  ## 2digit
str_split_fixed(str_split_fixed(x,":",2)[,2],"\\*",4)[,1:2] ## 4digit
