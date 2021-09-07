### HLA typing 결과 빈도 계산

setwd("~/Desktop/KCDC/HLAimputation/HLAtyping/all/")
ref <- read.table("~/Desktop/KCDC/HLAimputation/IMPUTE4/gen.calling.test/RESULTs/test/2digit.own.call.gen.fam")


head(ngs)
###
library(stringr)
ori <- read.csv("HLA.type.result.8genes.merged(2019.2020).csv")
head(ori)
colnames(ori) <- c("KID","YSample","A_1","A_2","B_1","B_2","C_1","C_2","DRB1_1","DRB1_2"
                   ,"DQA1_1","DQA1_2","DQB1_1","DQB1_2"
                   ,"DPA1_1","DPA1_2","DPB1_1","DPB1_2")
ngs <-ori
ngs$A.1 <- str_split_fixed(ngs$A_1,":",2)[,1]
ngs$A.1 <- str_split_fixed(ngs$A_1,":",2)[,1]
ngs$A.2 <- str_split_fixed(ngs$A_2,":",2)[,1]
ngs$B.1 <- str_split_fixed(ngs$B_1,":",2)[,1]
ngs$B.2 <- str_split_fixed(ngs$B_2,":",2)[,1]
ngs$C.1 <- str_split_fixed(ngs$C_1,":",2)[,1]
ngs$C.2 <- str_split_fixed(ngs$C_2,":",2)[,1]

ngs$DRB1.1 <- str_split_fixed(ngs$DRB1_1,":",2)[,1]
ngs$DRB1.2 <- str_split_fixed(ngs$DRB1_2,":",2)[,1]

ngs$DQA1.1 <- str_split_fixed(ngs$DQA1_1,":",2)[,1]
ngs$DQA1.2 <- str_split_fixed(ngs$DQA1_2,":",2)[,1]
ngs$DQB1.1 <- str_split_fixed(ngs$DQB1_1,":",2)[,1]
ngs$DQB1.2 <- str_split_fixed(ngs$DQB1_2,":",2)[,1]

ngs$DPA1.1 <- str_split_fixed(ngs$DPA1_1,":",2)[,1]
ngs$DPA1.2 <- str_split_fixed(ngs$DPA1_2,":",2)[,1]
ngs$DPB1.1 <- str_split_fixed(ngs$DPB1_1,":",2)[,1]
ngs$DPB1.2 <- str_split_fixed(ngs$DPB1_2,":",2)[,1]


colnames(ngs)


ngs.2d <- subset(ngs,select =c("KID","YSample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.2d)

###
ngs.2d[(!is.na(ngs.2d$DRB1.1)) & ((ngs.2d$DRB1.2 == ".")),]$DRB1.2 <- ngs.2d[(!is.na(ngs.2d$DRB1.1)) & (ngs.2d$DRB1.2 == "."),]$DRB1.1
ngs.2d[(ngs.2d$DRB1.1 == ".") & (!is.na(ngs.2d$DRB1.2)),]$DRB1.1 <- ngs.2d[(ngs.2d$DRB1.1== ".") & (!is.na(ngs.2d$DRB1.2)),]$DRB1.2
#ngs.2d[(!is.na(ngs.2d$DRB3.1)) & (is.na(ngs.2d$DRB3.2)),]$DRB3.2 <- ngs.2d[(!is.na(ngs.2d$DRB3.1)) & (is.na(ngs.2d$DRB3.2)),]$DRB3.1
#ngs.2d[(is.na(ngs.2d$DRB3.1)) & (!is.na(ngs.2d$DRB3.2)),]$DRB3.1 <- ngs.2d[(is.na(ngs.2d$DRB3.1)) & (!is.na(ngs.2d$DRB3.2)),]$DRB3.2

ngs.2d[(!is.na(ngs.2d$DQA1.1)) & (ngs.2d$DQA1.2== "."),]$DQA1.2 <- ngs.2d[(!is.na(ngs.2d$DQA1.1)) & (ngs.2d$DQA1.2== "."),]$DQA1.1
ngs.2d[(ngs.2d$DQA1.1== ".") & (!is.na(ngs.2d$DQA1.2)),]$DQA1.1 <- ngs.2d[(ngs.2d$DQA1.1== ".") & (!is.na(ngs.2d$DQA1.2)),]$DQA1.2
ngs.2d[(!is.na(ngs.2d$DQB1.1)) & (ngs.2d$DQB1.2== "."),]$DQB1.2 <- ngs.2d[(!is.na(ngs.2d$DQB1.1)) & (ngs.2d$DQB1.2== "."),]$DQB1.1
ngs.2d[(ngs.2d$DQB1.1== ".") & (!is.na(ngs.2d$DQB1.2)),]$DQB1.1 <- ngs.2d[(ngs.2d$DQB1.1== ".") & (!is.na(ngs.2d$DQB1.2)),]$DQB1.2

ngs.2d[(!is.na(ngs.2d$DPA1.1)) & (ngs.2d$DPA1.2== "."),]$DPA1.2 <- ngs.2d[(!is.na(ngs.2d$DPA1.1)) & (ngs.2d$DPA1.2== "."),]$DPA1.1
ngs.2d[(ngs.2d$DPA1.1== ".") & (!is.na(ngs.2d$DPA1.2)),]$DPA1.1 <- ngs.2d[(ngs.2d$DPA1.1== ".") & (!is.na(ngs.2d$DPA1.2)),]$DPA1.2
ngs.2d[(!is.na(ngs.2d$DPB1.1)) & (ngs.2d$DPB1.2== "."),]$DPB1.2 <- ngs.2d[(!is.na(ngs.2d$DPB1.1)) & (ngs.2d$DPB1.2== "."),]$DPB1.1
ngs.2d[(ngs.2d$DPB1.1== ".") & (!is.na(ngs.2d$DPB1.2)),]$DPB1.1 <- ngs.2d[(ngs.2d$DPB1.1== ".") & (!is.na(ngs.2d$DPB1.2)),]$DPB1.2

ngs.2d[(!is.na(ngs.2d$A.1)) & (ngs.2d$A.2== "."),]$A.2 <- ngs.2d[(!is.na(ngs.2d$A.1)) & (ngs.2d$A.2== "."),]$A.1
ngs.2d[(ngs.2d$A.1== ".") & (!is.na(ngs.2d$A.2)),]$A.1 <- ngs.2d[(ngs.2d$A.1== ".") & (!is.na(ngs.2d$A.2)),]$A.2

ngs.2d[(!is.na(ngs.2d$B.1)) & (ngs.2d$B.2== "."),]$B.2 <- ngs.2d[(!is.na(ngs.2d$B.1)) & (ngs.2d$B.2== "."),]$B.1
ngs.2d[(ngs.2d$B.1== ".") & (!is.na(ngs.2d$B.2)),]$B.1 <- ngs.2d[(ngs.2d$B.1== ".") & (!is.na(ngs.2d$B.2)),]$B.2

ngs.2d[(!is.na(ngs.2d$C.1)) & (ngs.2d$C.2== "."),]$C.2 <- ngs.2d[(!is.na(ngs.2d$C.1)) & (ngs.2d$C.2== "."),]$C.1
ngs.2d[(ngs.2d$C.1== ".") & (!is.na(ngs.2d$C.2)),]$C.1 <- ngs.2d[(ngs.2d$C.1== ".") & (!is.na(ngs.2d$C.2)),]$C.2


###





#a[(!is.na(a$DRB1)) & (is.na(a$DRB2)),]$DRB2 <- a[(!is.na(a$DRB1)) & (is.na(a$DRB2)),]$DRB1

head(ngs.2d)
ngs.2d

################4digit


colnames(ngs)
head(ngs)
##### 4digit
ngs$A.1 <- paste0(str_split_fixed(ngs$A_1,":",3)[,1],":",str_split_fixed(ngs$A_1,":",3)[,2])
ngs$A.2 <- paste0(str_split_fixed(ngs$A_2,":",3)[,1],":",str_split_fixed(ngs$A_2,":",3)[,2])
ngs$B.1 <- paste0(str_split_fixed(ngs$B_1,":",3)[,1],":",str_split_fixed(ngs$B_1,":",3)[,2])
ngs$B.2 <- paste0(str_split_fixed(ngs$B_2,":",3)[,1],":",str_split_fixed(ngs$B_2,":",3)[,2])
ngs$C.1 <- paste0(str_split_fixed(ngs$C_1,":",3)[,1],":",str_split_fixed(ngs$C_1,":",3)[,2])
ngs$C.2 <- paste0(str_split_fixed(ngs$C_2,":",3)[,1],":",str_split_fixed(ngs$C_2,":",3)[,2])

ngs$DRB1.1 <- paste0(str_split_fixed(ngs$DRB1_1,":",3)[,1],":",str_split_fixed(ngs$DRB1_1,":",3)[,2])
ngs$DRB1.2 <- paste0(str_split_fixed(ngs$DRB1_2,":",3)[,1],":",str_split_fixed(ngs$DRB1_2,":",3)[,2])

ngs$DQA1.1 <- paste0(str_split_fixed(ngs$DQA1_1,":",3)[,1],":",str_split_fixed(ngs$DQA1_1,":",3)[,2])
ngs$DQA1.2 <- paste0(str_split_fixed(ngs$DQA1_2,":",3)[,1],":",str_split_fixed(ngs$DQA1_2,":",3)[,2])
ngs$DQB1.1 <- paste0(str_split_fixed(ngs$DQB1_1,":",3)[,1],":",str_split_fixed(ngs$DQB1_1,":",3)[,2])
ngs$DQB1.2 <- paste0(str_split_fixed(ngs$DQB1_2,":",3)[,1],":",str_split_fixed(ngs$DQB1_2,":",3)[,2])

ngs$DPA1.1 <- paste0(str_split_fixed(ngs$DPA1_1,":",3)[,1],":",str_split_fixed(ngs$DPA1_1,":",3)[,2])
ngs$DPA1.2 <- paste0(str_split_fixed(ngs$DPA1_2,":",3)[,1],":",str_split_fixed(ngs$DPA1_2,":",3)[,2])
ngs$DPB1.1 <- paste0(str_split_fixed(ngs$DPB1_1,":",3)[,1],":",str_split_fixed(ngs$DPB1_1,":",3)[,2])
ngs$DPB1.2 <- paste0(str_split_fixed(ngs$DPB1_2,":",3)[,1],":",str_split_fixed(ngs$DPB1_2,":",3)[,2])



ngs.4d <- subset(ngs,select =c("KID","YSample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.4d)
head(ngs.4d)
ngs.4d[(!is.na(ngs.4d$DRB1.1)) & ((ngs.4d$DRB1.2 == ".:")),]$DRB1.2 <- ngs.4d[(!is.na(ngs.4d$DRB1.1)) & (ngs.4d$DRB1.2 == ".:"),]$DRB1.1
ngs.4d[(ngs.4d$DRB1.1 == ".:") & (!is.na(ngs.4d$DRB1.2)),]$DRB1.1 <- ngs.4d[(ngs.4d$DRB1.1== ".:") & (!is.na(ngs.4d$DRB1.2)),]$DRB1.2
#ngs.4d[(!is.na(ngs.4d$DRB3.1)) & (is.na(ngs.4d$DRB3.2)),]$DRB3.2 <- ngs.4d[(!is.na(ngs.4d$DRB3.1)) & (is.na(ngs.4d$DRB3.2)),]$DRB3.1
#ngs.4d[(is.na(ngs.4d$DRB3.1)) & (!is.na(ngs.4d$DRB3.2)),]$DRB3.1 <- ngs.4d[(is.na(ngs.4d$DRB3.1)) & (!is.na(ngs.4d$DRB3.2)),]$DRB3.2

ngs.4d[(!is.na(ngs.4d$DQA1.1)) & (ngs.4d$DQA1.2== ".:"),]$DQA1.2 <- ngs.4d[(!is.na(ngs.4d$DQA1.1)) & (ngs.4d$DQA1.2== ".:"),]$DQA1.1
ngs.4d[(ngs.4d$DQA1.1== ".:") & (!is.na(ngs.4d$DQA1.2)),]$DQA1.1 <- ngs.4d[(ngs.4d$DQA1.1== ".:") & (!is.na(ngs.4d$DQA1.2)),]$DQA1.2
ngs.4d[(!is.na(ngs.4d$DQB1.1)) & (ngs.4d$DQB1.2== ".:"),]$DQB1.2 <- ngs.4d[(!is.na(ngs.4d$DQB1.1)) & (ngs.4d$DQB1.2== ".:"),]$DQB1.1
ngs.4d[(ngs.4d$DQB1.1== ".:") & (!is.na(ngs.4d$DQB1.2)),]$DQB1.1 <- ngs.4d[(ngs.4d$DQB1.1== ".:") & (!is.na(ngs.4d$DQB1.2)),]$DQB1.2

ngs.4d[(!is.na(ngs.4d$DPA1.1)) & (ngs.4d$DPA1.2== ".:"),]$DPA1.2 <- ngs.4d[(!is.na(ngs.4d$DPA1.1)) & (ngs.4d$DPA1.2== ".:"),]$DPA1.1
ngs.4d[(ngs.4d$DPA1.1== ".:") & (!is.na(ngs.4d$DPA1.2)),]$DPA1.1 <- ngs.4d[(ngs.4d$DPA1.1== ".:") & (!is.na(ngs.4d$DPA1.2)),]$DPA1.2
ngs.4d[(!is.na(ngs.4d$DPB1.1)) & (ngs.4d$DPB1.2== ".:"),]$DPB1.2 <- ngs.4d[(!is.na(ngs.4d$DPB1.1)) & (ngs.4d$DPB1.2== ".:"),]$DPB1.1
ngs.4d[(ngs.4d$DPB1.1== ".:") & (!is.na(ngs.4d$DPB1.2)),]$DPB1.1 <- ngs.4d[(ngs.4d$DPB1.1== ".:") & (!is.na(ngs.4d$DPB1.2)),]$DPB1.2

ngs.4d[(!is.na(ngs.4d$A.1)) & (ngs.4d$A.2== ".:"),]$A.2 <- ngs.4d[(!is.na(ngs.4d$A.1)) & (ngs.4d$A.2== ".:"),]$A.1
ngs.4d[(ngs.4d$A.1== ".:") & (!is.na(ngs.4d$A.2)),]$A.1 <- ngs.4d[(ngs.4d$A.1== ".:") & (!is.na(ngs.4d$A.2)),]$A.2

ngs.4d[(!is.na(ngs.4d$B.1)) & (ngs.4d$B.2== ".:"),]$B.2 <- ngs.4d[(!is.na(ngs.4d$B.1)) & (ngs.4d$B.2== ".:"),]$B.1
ngs.4d[(ngs.4d$B.1== ".:") & (!is.na(ngs.4d$B.2)),]$B.1 <- ngs.4d[(ngs.4d$B.1== ".:") & (!is.na(ngs.4d$B.2)),]$B.2

ngs.4d[(!is.na(ngs.4d$C.1)) & (ngs.4d$C.2== ".:"),]$C.2 <- ngs.4d[(!is.na(ngs.4d$C.1)) & (ngs.4d$C.2== ".:"),]$C.1
ngs.4d[(ngs.4d$C.1== ".:") & (!is.na(ngs.4d$C.2)),]$C.1 <- ngs.4d[(ngs.4d$C.1== ".:") & (!is.na(ngs.4d$C.2)),]$C.2


head(ngs.4d)





head(ngs.2d)
head(ngs.4d)
ngs.2d <- ngs.2d[ngs.2d$KID %in% ref$V1,]
ngs.4d <- ngs.4d[ngs.4d$KID %in% ref$V1,]



for (i in 3:18) {
  ngs.2d[,i] <- str_replace_all(ngs.2d[,i],"\\.","NULL_allele")
  ngs.2d[,i] <- str_replace_all(ngs.2d[,i],"","NULL_allele")
  #ngs.2d[,i] <- str_replace_na(ngs.2d[,i],"")
  
  #ngs.2d[,i] <- ngs.2d[,i] %>% replace_na('')
  #  df[,i] <- str_replace_all(df[,i],is.na(),"NULL_allele")
}


for (i in 3:18) {
  ngs.4d[,i] <- str_replace_all(ngs.4d[,i],"\\.:","NULL_allele")
  ngs.4d[,i] <- str_replace_all(ngs.4d[,i],"\\.","NULL_allele")
  ngs.4d[,i] <- str_replace_all(ngs.4d[,i],"",":")
  #ngs.4d[,i] <- str_replace_na(ngs.4d[,i],"")
  #  df[,i] <- str_replace_all(df[,i],is.na(),"NULL_allele")
}



first_allele <- ngs.2d[,c("KID","YSample","A.1","B.1","C.1","DRB1.1","DQA1.1","DQB1.1","DPA1.1","DPB1.1")]
second_allele <- ngs.2d[,c("KID","YSample","A.2","B.2","C.2","DRB1.2","DQA1.2","DQB1.2","DPA1.2","DPB1.2")]
head(first_allele)
head(second_allele)
colnames(first_allele) <- c("KID","YSample","A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1")
colnames(second_allele) <- c("KID","YSample","A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1")
out <- rbind(first_allele,second_allele)


A <- as.data.frame(table(out$A))
colnames(A)<-c("HLA_type","Freq")
A$HLA_gene <- "HLA-A"

B <- as.data.frame(table(out$B))
colnames(B)<-c("HLA_type","Freq")
B$HLA_gene <- "HLA-B"

C <- as.data.frame(table(out$C))
colnames(C)<-c("HLA_type","Freq")
C$HLA_gene <- "HLA-C"

DRB1 <- as.data.frame(table(out$DRB1))
colnames(DRB1)<-c("HLA_type","Freq")
DRB1$HLA_gene <- "HLA-DRB1"

DQA1 <- as.data.frame(table(out$DQA1))
colnames(DQA1)<-c("HLA_type","Freq")
DQA1$HLA_gene <- "HLA-DQA1"

DQB1 <- as.data.frame(table(out$DQB1))
colnames(DQB1)<-c("HLA_type","Freq")
DQB1$HLA_gene <- "HLA-DQB1"

DPA1 <- as.data.frame(table(out$DPA1))
colnames(DPA1)<-c("HLA_type","Freq")
DPA1$HLA_gene <- "HLA-DPA1"

DPB1 <- as.data.frame(table(out$DPB1))
colnames(DPB1)<-c("HLA_type","Freq")
DPB1$HLA_gene <- "HLA-DPB1"

td <- rbind(A,B)
td <- rbind(td,C)
td <- rbind(td,DRB1)
td <- rbind(td,DQA1)
td <- rbind(td,DQB1)
td <- rbind(td,DPA1)
td <- rbind(td,DPB1)
head(td)
td$digit <- '2digit'
td$freq <- td$Freq/1022
td <- td[,c('digit','HLA_gene','HLA_type','freq')]


first_allele <- ngs.4d[,c("KID","YSample","A.1","B.1","C.1","DRB1.1","DQA1.1","DQB1.1","DPA1.1","DPB1.1")]
second_allele <- ngs.4d[,c("KID","YSample","A.2","B.2","C.2","DRB1.2","DQA1.2","DQB1.2","DPA1.2","DPB1.2")]
head(first_allele)
head(second_allele)
colnames(first_allele) <- c("KID","YSample","A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1")
colnames(second_allele) <- c("KID","YSample","A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1")
out <- rbind(first_allele,second_allele)


A <- as.data.frame(table(out$A))
colnames(A)<-c("HLA_type","Freq")
A$HLA_gene <- "HLA-A"

B <- as.data.frame(table(out$B))
colnames(B)<-c("HLA_type","Freq")
B$HLA_gene <- "HLA-B"

C <- as.data.frame(table(out$C))
colnames(C)<-c("HLA_type","Freq")
C$HLA_gene <- "HLA-C"

DRB1 <- as.data.frame(table(out$DRB1))
colnames(DRB1)<-c("HLA_type","Freq")
DRB1$HLA_gene <- "HLA-DRB1"

DQA1 <- as.data.frame(table(out$DQA1))
colnames(DQA1)<-c("HLA_type","Freq")
DQA1$HLA_gene <- "HLA-DQA1"

DQB1 <- as.data.frame(table(out$DQB1))
colnames(DQB1)<-c("HLA_type","Freq")
DQB1$HLA_gene <- "HLA-DQB1"

DPA1 <- as.data.frame(table(out$DPA1))
colnames(DPA1)<-c("HLA_type","Freq")
DPA1$HLA_gene <- "HLA-DPA1"

DPB1 <- as.data.frame(table(out$DPB1))
colnames(DPB1)<-c("HLA_type","Freq")
DPB1$HLA_gene <- "HLA-DPB1"

fd <- rbind(A,B)
fd <- rbind(fd,C)
fd <- rbind(fd,DRB1)
fd <- rbind(fd,DQA1)
fd <- rbind(fd,DQB1)
fd <- rbind(fd,DPA1)
fd <- rbind(fd,DPB1)
head(fd)
fd$digit <- '4digit'
fd$freq <- fd$Freq/1022
fd <- fd[,c('digit','HLA_gene','HLA_type','freq')]


head(td)
head(fd)


df <- rbind(td,fd)
head(df)
library(tidyr)
#df$HLA_type <- df$HLA_type %>% replace_na('NULL_allele')
head(df)

write.table(df,"HLA_type_frequency.txt",col.names = T,row.names = F,quote = F,na = "NULL_allele",sep = "\t")


#박사님 HLA typing 결과 frequency 계산값입니다.
#511 샘플에 대해서 진행했습니다.(1022 allele을 나누어 진행)
#NULL allele에 대한 frequency 도 포함하였습니다.