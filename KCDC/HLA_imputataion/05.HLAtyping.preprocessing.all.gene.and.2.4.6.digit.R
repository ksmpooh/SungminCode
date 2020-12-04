###########################20200828
###HLA typing 결과 정리
#2019
setwd("c:/Users/user/Desktop/KCDC/transplantation/HLAtyping/")
ori <- read.csv("HLAtyping.allGene.result.csv",header = T)
#2020
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/HLAtyping/2020/")
ori <- read.csv("2020_HLAtyping_tocheck.csv",header = T)
tail(ori)
head(ori)
colnames(ori) <- c("X","A_1","A_2","B_1","B_2","C_1","C_2","DRB1_1","DRB1_2"
                   ,"DRB3_1","DRB3_2","DRB4_1","DRB4_2","DRB5_1","DRB5_2"
                   ,"DQA1_1","DQA1_2","DQB1_1","DQB1_2"
                   ,"DPA1_1","DPA1_2","DPB1_1","DPB1_2")
####
library(stringr)





colnames(ori)
head(ori)

ngs <-ori
head(ngs)
head(ngs)
###############2 digit
ngs$A.1 <- str_split_fixed(str_split_fixed(ngs$A_1,":",2)[,1],"\\*",2)[,2]
ngs$A.2 <- str_split_fixed(str_split_fixed(ngs$A_2,":",2)[,1],"\\*",2)[,2]
ngs$B.1 <- str_split_fixed(str_split_fixed(ngs$B_1,":",2)[,1],"\\*",2)[,2]
ngs$B.2 <- str_split_fixed(str_split_fixed(ngs$B_2,":",2)[,1],"\\*",2)[,2]
ngs$C.1 <- str_split_fixed(str_split_fixed(ngs$C_1,":",2)[,1],"\\*",2)[,2]
ngs$C.2 <- str_split_fixed(str_split_fixed(ngs$C_2,":",2)[,1],"\\*",2)[,2]

ngs$DRB1.1 <- str_split_fixed(str_split_fixed(ngs$DRB1_1,":",2)[,1],"\\*",2)[,2]
ngs$DRB1.2 <- str_split_fixed(str_split_fixed(ngs$DRB1_2,":",2)[,1],"\\*",2)[,2]
ngs$DRB3.1 <- str_split_fixed(str_split_fixed(ngs$DRB3_1,":",2)[,1],"\\*",2)[,2]
ngs$DRB3.2 <- str_split_fixed(str_split_fixed(ngs$DRB3_2,":",2)[,1],"\\*",2)[,2]

ngs$DQA1.1 <- str_split_fixed(str_split_fixed(ngs$DQA1_1,":",2)[,1],"\\*",2)[,2]
ngs$DQA1.2 <- str_split_fixed(str_split_fixed(ngs$DQA1_2,":",2)[,1],"\\*",2)[,2]
ngs$DQB1.1 <- str_split_fixed(str_split_fixed(ngs$DQB1_1,":",2)[,1],"\\*",2)[,2]
ngs$DQB1.2 <- str_split_fixed(str_split_fixed(ngs$DQB1_2,":",2)[,1],"\\*",2)[,2]

ngs$DPA1.1 <- str_split_fixed(str_split_fixed(ngs$DPA1_1,":",2)[,1],"\\*",2)[,2]
ngs$DPA1.2 <- str_split_fixed(str_split_fixed(ngs$DPA1_2,":",2)[,1],"\\*",2)[,2]
ngs$DPB1.1 <- str_split_fixed(str_split_fixed(ngs$DPB1_1,":",2)[,1],"\\*",2)[,2]
ngs$DPB1.2 <- str_split_fixed(str_split_fixed(ngs$DPB1_2,":",2)[,1],"\\*",2)[,2]


colnames(ngs)


ngs.2d <- subset(ngs,select =c("X","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DRB3.1","DRB3.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.2d)
features <- colnames(ngs.2d)[2:ncol(ngs.2d)]
features
for (i in features){
  #ngs.2d[,i] <- strtoi(ngs.2d[,i])
  ngs.2d[,i] <- as.integer(ngs.2d[,i])
}

ngs.2d[(!is.na(ngs.2d$DRB1.1)) & (is.na(ngs.2d$DRB1.2)),]$DRB1.2 <- ngs.2d[(!is.na(ngs.2d$DRB1.1)) & (is.na(ngs.2d$DRB1.2)),]$DRB1.1
ngs.2d[(is.na(ngs.2d$DRB1.1)) & (!is.na(ngs.2d$DRB1.2)),]$DRB1.1 <- ngs.2d[(is.na(ngs.2d$DRB1.1)) & (!is.na(ngs.2d$DRB1.2)),]$DRB1.2
ngs.2d[(!is.na(ngs.2d$DRB3.1)) & (is.na(ngs.2d$DRB3.2)),]$DRB3.2 <- ngs.2d[(!is.na(ngs.2d$DRB3.1)) & (is.na(ngs.2d$DRB3.2)),]$DRB3.1
ngs.2d[(is.na(ngs.2d$DRB3.1)) & (!is.na(ngs.2d$DRB3.2)),]$DRB3.1 <- ngs.2d[(is.na(ngs.2d$DRB3.1)) & (!is.na(ngs.2d$DRB3.2)),]$DRB3.2

ngs.2d[(!is.na(ngs.2d$DQA1.1)) & (is.na(ngs.2d$DQA1.2)),]$DQA1.2 <- ngs.2d[(!is.na(ngs.2d$DQA1.1)) & (is.na(ngs.2d$DQA1.2)),]$DQA1.1
ngs.2d[(is.na(ngs.2d$DQA1.1)) & (!is.na(ngs.2d$DQA1.2)),]$DQA1.1 <- ngs.2d[(is.na(ngs.2d$DQA1.1)) & (!is.na(ngs.2d$DQA1.2)),]$DQA1.2
ngs.2d[(!is.na(ngs.2d$DQB1.1)) & (is.na(ngs.2d$DQB1.2)),]$DQB1.2 <- ngs.2d[(!is.na(ngs.2d$DQB1.1)) & (is.na(ngs.2d$DQB1.2)),]$DQB1.1
ngs.2d[(is.na(ngs.2d$DQB1.1)) & (!is.na(ngs.2d$DQB1.2)),]$DQB1.1 <- ngs.2d[(is.na(ngs.2d$DQB1.1)) & (!is.na(ngs.2d$DQB1.2)),]$DQB1.2

ngs.2d[(!is.na(ngs.2d$DPA1.1)) & (is.na(ngs.2d$DPA1.2)),]$DPA1.2 <- ngs.2d[(!is.na(ngs.2d$DPA1.1)) & (is.na(ngs.2d$DPA1.2)),]$DPA1.1
ngs.2d[(is.na(ngs.2d$DPA1.1)) & (!is.na(ngs.2d$DPA1.2)),]$DPA1.1 <- ngs.2d[(is.na(ngs.2d$DPA1.1)) & (!is.na(ngs.2d$DPA1.2)),]$DPA1.2
ngs.2d[(!is.na(ngs.2d$DPB1.1)) & (is.na(ngs.2d$DPB1.2)),]$DPB1.2 <- ngs.2d[(!is.na(ngs.2d$DPB1.1)) & (is.na(ngs.2d$DPB1.2)),]$DPB1.1
ngs.2d[(is.na(ngs.2d$DPB1.1)) & (!is.na(ngs.2d$DPB1.2)),]$DPB1.1 <- ngs.2d[(is.na(ngs.2d$DPB1.1)) & (!is.na(ngs.2d$DPB1.2)),]$DPB1.2

ngs.2d[(!is.na(ngs.2d$A.1)) & (is.na(ngs.2d$A.2)),]$A.2 <- ngs.2d[(!is.na(ngs.2d$A.1)) & (is.na(ngs.2d$A.2)),]$A.1
ngs.2d[(is.na(ngs.2d$A.1)) & (!is.na(ngs.2d$A.2)),]$A.1 <- ngs.2d[(is.na(ngs.2d$A.1)) & (!is.na(ngs.2d$A.2)),]$A.2

ngs.2d[(!is.na(ngs.2d$B.1)) & (is.na(ngs.2d$B.2)),]$B.2 <- ngs.2d[(!is.na(ngs.2d$B.1)) & (is.na(ngs.2d$B.2)),]$B.1
ngs.2d[(is.na(ngs.2d$B.1)) & (!is.na(ngs.2d$B.2)),]$B.1 <- ngs.2d[(is.na(ngs.2d$B.1)) & (!is.na(ngs.2d$B.2)),]$B.2

ngs.2d[(!is.na(ngs.2d$C.1)) & (is.na(ngs.2d$C.2)),]$C.2 <- ngs.2d[(!is.na(ngs.2d$C.1)) & (is.na(ngs.2d$C.2)),]$C.1
ngs.2d[(is.na(ngs.2d$C.1)) & (!is.na(ngs.2d$C.2)),]$C.1 <- ngs.2d[(is.na(ngs.2d$C.1)) & (!is.na(ngs.2d$C.2)),]$C.2

#a[(!is.na(a$DRB1)) & (is.na(a$DRB2)),]$DRB2 <- a[(!is.na(a$DRB1)) & (is.na(a$DRB2)),]$DRB1



################4digit


colnames(ngs)
head(ngs)
##### 4digit
ngs$A.1 <- str_split_fixed(paste0(str_split_fixed(ngs$A_1,":",3)[,1],str_split_fixed(ngs$A_1,":",3)[,2]),"\\*",2)[,2]
ngs$A.2 <- str_split_fixed(paste0(str_split_fixed(ngs$A_2,":",3)[,1],str_split_fixed(ngs$A_2,":",3)[,2]),"\\*",2)[,2]
ngs$B.1 <- str_split_fixed(paste0(str_split_fixed(ngs$B_1,":",3)[,1],str_split_fixed(ngs$B_1,":",3)[,2]),"\\*",2)[,2]
ngs$B.2 <- str_split_fixed(paste0(str_split_fixed(ngs$B_2,":",3)[,1],str_split_fixed(ngs$B_2,":",3)[,2]),"\\*",2)[,2]
ngs$C.1 <- str_split_fixed(paste0(str_split_fixed(ngs$C_1,":",3)[,1],str_split_fixed(ngs$C_1,":",3)[,2]),"\\*",2)[,2]
ngs$C.2 <- str_split_fixed(paste0(str_split_fixed(ngs$C_2,":",3)[,1],str_split_fixed(ngs$C_2,":",3)[,2]),"\\*",2)[,2]

ngs$DRB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB1_1,":",3)[,1],str_split_fixed(ngs$DRB1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DRB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB1_2,":",3)[,1],str_split_fixed(ngs$DRB1_2,":",3)[,2]),"\\*",2)[,2]
ngs$DRB3.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB3_1,":",3)[,1],str_split_fixed(ngs$DRB3_1,":",3)[,2]),"\\*",2)[,2]
ngs$DRB3.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB3_2,":",3)[,1],str_split_fixed(ngs$DRB3_2,":",3)[,2]),"\\*",2)[,2]

ngs$DQA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_1,":",3)[,1],str_split_fixed(ngs$DQA1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DQA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_2,":",3)[,1],str_split_fixed(ngs$DQA1_2,":",3)[,2]),"\\*",2)[,2]
ngs$DQB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_1,":",3)[,1],str_split_fixed(ngs$DQB1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DQB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_2,":",3)[,1],str_split_fixed(ngs$DQB1_2,":",3)[,2]),"\\*",2)[,2]

ngs$DPA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_1,":",3)[,1],str_split_fixed(ngs$DPA1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DPA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_2,":",3)[,1],str_split_fixed(ngs$DPA1_2,":",3)[,2]),"\\*",2)[,2]
ngs$DPB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_1,":",3)[,1],str_split_fixed(ngs$DPB1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DPB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_2,":",3)[,1],str_split_fixed(ngs$DPB1_2,":",3)[,2]),"\\*",2)[,2]



ngs.4d <- subset(ngs,select =c("X","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DRB3.1","DRB3.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.4d)
features <- colnames(ngs.4d)[2:ncol(ngs.4d)]
features
for (i in features){
  #ngs.4d[,i] <- strtoi(ngs.4d[,i])
  ngs.4d[,i] <- as.integer(ngs.4d[,i])
}

ngs.4d[(!is.na(ngs.4d$DRB1.1)) & (is.na(ngs.4d$DRB1.2)),]$DRB1.2 <- ngs.4d[(!is.na(ngs.4d$DRB1.1)) & (is.na(ngs.4d$DRB1.2)),]$DRB1.1
ngs.4d[(is.na(ngs.4d$DRB1.1)) & (!is.na(ngs.4d$DRB1.2)),]$DRB1.1 <- ngs.4d[(is.na(ngs.4d$DRB1.1)) & (!is.na(ngs.4d$DRB1.2)),]$DRB1.2
ngs.4d[(!is.na(ngs.4d$DRB3.1)) & (is.na(ngs.4d$DRB3.2)),]$DRB3.2 <- ngs.4d[(!is.na(ngs.4d$DRB3.1)) & (is.na(ngs.4d$DRB3.2)),]$DRB3.1
ngs.4d[(is.na(ngs.4d$DRB3.1)) & (!is.na(ngs.4d$DRB3.2)),]$DRB3.1 <- ngs.4d[(is.na(ngs.4d$DRB3.1)) & (!is.na(ngs.4d$DRB3.2)),]$DRB3.2

ngs.4d[(!is.na(ngs.4d$DQA1.1)) & (is.na(ngs.4d$DQA1.2)),]$DQA1.2 <- ngs.4d[(!is.na(ngs.4d$DQA1.1)) & (is.na(ngs.4d$DQA1.2)),]$DQA1.1
ngs.4d[(is.na(ngs.4d$DQA1.1)) & (!is.na(ngs.4d$DQA1.2)),]$DQA1.1 <- ngs.4d[(is.na(ngs.4d$DQA1.1)) & (!is.na(ngs.4d$DQA1.2)),]$DQA1.2
ngs.4d[(!is.na(ngs.4d$DQB1.1)) & (is.na(ngs.4d$DQB1.2)),]$DQB1.2 <- ngs.4d[(!is.na(ngs.4d$DQB1.1)) & (is.na(ngs.4d$DQB1.2)),]$DQB1.1
ngs.4d[(is.na(ngs.4d$DQB1.1)) & (!is.na(ngs.4d$DQB1.2)),]$DQB1.1 <- ngs.4d[(is.na(ngs.4d$DQB1.1)) & (!is.na(ngs.4d$DQB1.2)),]$DQB1.2

ngs.4d[(!is.na(ngs.4d$DPA1.1)) & (is.na(ngs.4d$DPA1.2)),]$DPA1.2 <- ngs.4d[(!is.na(ngs.4d$DPA1.1)) & (is.na(ngs.4d$DPA1.2)),]$DPA1.1
ngs.4d[(is.na(ngs.4d$DPA1.1)) & (!is.na(ngs.4d$DPA1.2)),]$DPA1.1 <- ngs.4d[(is.na(ngs.4d$DPA1.1)) & (!is.na(ngs.4d$DPA1.2)),]$DPA1.2
ngs.4d[(!is.na(ngs.4d$DPB1.1)) & (is.na(ngs.4d$DPB1.2)),]$DPB1.2 <- ngs.4d[(!is.na(ngs.4d$DPB1.1)) & (is.na(ngs.4d$DPB1.2)),]$DPB1.1
ngs.4d[(is.na(ngs.4d$DPB1.1)) & (!is.na(ngs.4d$DPB1.2)),]$DPB1.1 <- ngs.4d[(is.na(ngs.4d$DPB1.1)) & (!is.na(ngs.4d$DPB1.2)),]$DPB1.2

ngs.4d[(!is.na(ngs.4d$A.1)) & (is.na(ngs.4d$A.2)),]$A.2 <- ngs.4d[(!is.na(ngs.4d$A.1)) & (is.na(ngs.4d$A.2)),]$A.1
ngs.4d[(is.na(ngs.4d$A.1)) & (!is.na(ngs.4d$A.2)),]$A.1 <- ngs.4d[(is.na(ngs.4d$A.1)) & (!is.na(ngs.4d$A.2)),]$A.2

ngs.4d[(!is.na(ngs.4d$B.1)) & (is.na(ngs.4d$B.2)),]$B.2 <- ngs.4d[(!is.na(ngs.4d$B.1)) & (is.na(ngs.4d$B.2)),]$B.1
ngs.4d[(is.na(ngs.4d$B.1)) & (!is.na(ngs.4d$B.2)),]$B.1 <- ngs.4d[(is.na(ngs.4d$B.1)) & (!is.na(ngs.4d$B.2)),]$B.2

ngs.4d[(!is.na(ngs.4d$C.1)) & (is.na(ngs.4d$C.2)),]$C.2 <- ngs.4d[(!is.na(ngs.4d$C.1)) & (is.na(ngs.4d$C.2)),]$C.1
ngs.4d[(is.na(ngs.4d$C.1)) & (!is.na(ngs.4d$C.2)),]$C.1 <- ngs.4d[(is.na(ngs.4d$C.1)) & (!is.na(ngs.4d$C.2)),]$C.2

head(ngs.4d)

######### 6digit
head(ngs)

ngs$A.1 <- str_split_fixed(paste0(str_split_fixed(ngs$A_1,":",4)[,1],str_split_fixed(ngs$A_1,":",4)[,2],str_split_fixed(ngs$A_1,":",4)[,3]),"\\*",2)[,2]
ngs$A.2 <- str_split_fixed(paste0(str_split_fixed(ngs$A_2,":",4)[,1],str_split_fixed(ngs$A_2,":",4)[,2],str_split_fixed(ngs$A_2,":",4)[,3]),"\\*",2)[,2]

ngs$B.1 <- str_split_fixed(paste0(str_split_fixed(ngs$B_1,":",4)[,1],str_split_fixed(ngs$B_1,":",4)[,2],str_split_fixed(ngs$B_1,":",4)[,3]),"\\*",2)[,2]
ngs$B.2 <- str_split_fixed(paste0(str_split_fixed(ngs$B_2,":",4)[,1],str_split_fixed(ngs$B_2,":",4)[,2],str_split_fixed(ngs$B_2,":",4)[,3]),"\\*",2)[,2]

ngs$C.1 <- str_split_fixed(paste0(str_split_fixed(ngs$C_1,":",4)[,1],str_split_fixed(ngs$C_1,":",4)[,2],str_split_fixed(ngs$C_1,":",4)[,3]),"\\*",2)[,2]
ngs$C.2 <- str_split_fixed(paste0(str_split_fixed(ngs$C_2,":",4)[,1],str_split_fixed(ngs$C_2,":",4)[,2],str_split_fixed(ngs$C_2,":",4)[,3]),"\\*",2)[,2]

ngs$DRB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB1_1,":",4)[,1],str_split_fixed(ngs$DRB1_1,":",4)[,2],str_split_fixed(ngs$DRB1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DRB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB1_2,":",4)[,1],str_split_fixed(ngs$DRB1_2,":",4)[,2],str_split_fixed(ngs$DRB1_2,":",4)[,3]),"\\*",2)[,2]
ngs$DRB3.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB3_1,":",4)[,1],str_split_fixed(ngs$DRB3_1,":",4)[,2],str_split_fixed(ngs$DRB3_1,":",4)[,3]),"\\*",2)[,2]
ngs$DRB3.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DRB3_2,":",4)[,1],str_split_fixed(ngs$DRB3_2,":",4)[,2],str_split_fixed(ngs$DRB3_2,":",4)[,3]),"\\*",2)[,2]

ngs$DQB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_1,":",4)[,1],str_split_fixed(ngs$DQB1_1,":",4)[,2],str_split_fixed(ngs$DQB1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DQB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_2,":",4)[,1],str_split_fixed(ngs$DQB1_2,":",4)[,2],str_split_fixed(ngs$DQB1_2,":",4)[,3]),"\\*",2)[,2]
ngs$DQA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_1,":",4)[,1],str_split_fixed(ngs$DQA1_1,":",4)[,2],str_split_fixed(ngs$DQA1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DQA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_2,":",4)[,1],str_split_fixed(ngs$DQA1_2,":",4)[,2],str_split_fixed(ngs$DQA1_2,":",4)[,3]),"\\*",2)[,2]

ngs$DPA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_1,":",4)[,1],str_split_fixed(ngs$DPA1_1,":",4)[,2],str_split_fixed(ngs$DPA1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DPA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_2,":",4)[,1],str_split_fixed(ngs$DPA1_2,":",4)[,2],str_split_fixed(ngs$DPA1_2,":",4)[,3]),"\\*",2)[,2]
ngs$DPB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_1,":",4)[,1],str_split_fixed(ngs$DPB1_1,":",4)[,2],str_split_fixed(ngs$DPB1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DPB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_2,":",4)[,1],str_split_fixed(ngs$DPB1_2,":",4)[,2],str_split_fixed(ngs$DPB1_2,":",4)[,3]),"\\*",2)[,2]

ngs.6d <- subset(ngs,select =c("X","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DRB3.1","DRB3.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.6d)
features <- colnames(ngs.6d)[2:ncol(ngs.6d)]
features
for (i in features){
  #ngs.6d[,i] <- strtoi(ngs.6d[,i])
  ngs.6d[,i] <- as.integer(ngs.6d[,i])
}

ngs.6d[(!is.na(ngs.6d$DRB1.1)) & (is.na(ngs.6d$DRB1.2)),]$DRB1.2 <- ngs.6d[(!is.na(ngs.6d$DRB1.1)) & (is.na(ngs.6d$DRB1.2)),]$DRB1.1
ngs.6d[(is.na(ngs.6d$DRB1.1)) & (!is.na(ngs.6d$DRB1.2)),]$DRB1.1 <- ngs.6d[(is.na(ngs.6d$DRB1.1)) & (!is.na(ngs.6d$DRB1.2)),]$DRB1.2
ngs.6d[(!is.na(ngs.6d$DRB3.1)) & (is.na(ngs.6d$DRB3.2)),]$DRB3.2 <- ngs.6d[(!is.na(ngs.6d$DRB3.1)) & (is.na(ngs.6d$DRB3.2)),]$DRB3.1
ngs.6d[(is.na(ngs.6d$DRB3.1)) & (!is.na(ngs.6d$DRB3.2)),]$DRB3.1 <- ngs.6d[(is.na(ngs.6d$DRB3.1)) & (!is.na(ngs.6d$DRB3.2)),]$DRB3.2

ngs.6d[(!is.na(ngs.6d$DQA1.1)) & (is.na(ngs.6d$DQA1.2)),]$DQA1.2 <- ngs.6d[(!is.na(ngs.6d$DQA1.1)) & (is.na(ngs.6d$DQA1.2)),]$DQA1.1
ngs.6d[(is.na(ngs.6d$DQA1.1)) & (!is.na(ngs.6d$DQA1.2)),]$DQA1.1 <- ngs.6d[(is.na(ngs.6d$DQA1.1)) & (!is.na(ngs.6d$DQA1.2)),]$DQA1.2
ngs.6d[(!is.na(ngs.6d$DQB1.1)) & (is.na(ngs.6d$DQB1.2)),]$DQB1.2 <- ngs.6d[(!is.na(ngs.6d$DQB1.1)) & (is.na(ngs.6d$DQB1.2)),]$DQB1.1
ngs.6d[(is.na(ngs.6d$DQB1.1)) & (!is.na(ngs.6d$DQB1.2)),]$DQB1.1 <- ngs.6d[(is.na(ngs.6d$DQB1.1)) & (!is.na(ngs.6d$DQB1.2)),]$DQB1.2

ngs.6d[(!is.na(ngs.6d$DPA1.1)) & (is.na(ngs.6d$DPA1.2)),]$DPA1.2 <- ngs.6d[(!is.na(ngs.6d$DPA1.1)) & (is.na(ngs.6d$DPA1.2)),]$DPA1.1
ngs.6d[(is.na(ngs.6d$DPA1.1)) & (!is.na(ngs.6d$DPA1.2)),]$DPA1.1 <- ngs.6d[(is.na(ngs.6d$DPA1.1)) & (!is.na(ngs.6d$DPA1.2)),]$DPA1.2
ngs.6d[(!is.na(ngs.6d$DPB1.1)) & (is.na(ngs.6d$DPB1.2)),]$DPB1.2 <- ngs.6d[(!is.na(ngs.6d$DPB1.1)) & (is.na(ngs.6d$DPB1.2)),]$DPB1.1
ngs.6d[(is.na(ngs.6d$DPB1.1)) & (!is.na(ngs.6d$DPB1.2)),]$DPB1.1 <- ngs.6d[(is.na(ngs.6d$DPB1.1)) & (!is.na(ngs.6d$DPB1.2)),]$DPB1.2

ngs.6d[(!is.na(ngs.6d$A.1)) & (is.na(ngs.6d$A.2)),]$A.2 <- ngs.6d[(!is.na(ngs.6d$A.1)) & (is.na(ngs.6d$A.2)),]$A.1
ngs.6d[(is.na(ngs.6d$A.1)) & (!is.na(ngs.6d$A.2)),]$A.1 <- ngs.6d[(is.na(ngs.6d$A.1)) & (!is.na(ngs.6d$A.2)),]$A.2

ngs.6d[(!is.na(ngs.6d$B.1)) & (is.na(ngs.6d$B.2)),]$B.2 <- ngs.6d[(!is.na(ngs.6d$B.1)) & (is.na(ngs.6d$B.2)),]$B.1
ngs.6d[(is.na(ngs.6d$B.1)) & (!is.na(ngs.6d$B.2)),]$B.1 <- ngs.6d[(is.na(ngs.6d$B.1)) & (!is.na(ngs.6d$B.2)),]$B.2

ngs.6d[(!is.na(ngs.6d$C.1)) & (is.na(ngs.6d$C.2)),]$C.2 <- ngs.6d[(!is.na(ngs.6d$C.1)) & (is.na(ngs.6d$C.2)),]$C.1
ngs.6d[(is.na(ngs.6d$C.1)) & (!is.na(ngs.6d$C.2)),]$C.1 <- ngs.6d[(is.na(ngs.6d$C.1)) & (!is.na(ngs.6d$C.2)),]$C.2

head(ngs.6d)
#######################################################

#2019
write.csv(ngs.2d,"20200828/HLAtyping.alle.gene.2digit.csv",row.names = F,quote = F)
write.csv(ngs.4d,"20200828/HLAtyping.alle.gene.4digit.csv",row.names = F,quote = F)
write.csv(ngs.6d,"20200828/HLAtyping.alle.gene.6digit.csv",row.names = F,quote = F)


#2020
write.csv(ngs.2d,"HLAtyping.alle.gene.2digit.csv",row.names = F,quote = F)
write.csv(ngs.4d,"HLAtyping.alle.gene.4digit.csv",row.names = F,quote = F)



####################add kchip ID
##2019
setwd("c:/Users/user/Desktop/KCDC/transplantation/HLAtyping/")
ref <- read.csv("HLA_NGS_typing_255samples_results_202002.csv")
head(ref)
df <- read.csv("20200828/HLAtyping.alle.gene.2digit.csv")
df <- merge(ref[,c(1,2)],df,by.x = "YSample",by.y = "X")
colnames(df)
colnames(df) <- c("YSample","KID","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2","NGS_DRB3.1","NGS_DRB3.2","NGS_DQA1.1", 
                   "NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
write.csv(df,"20200828/HLAtyping.alle.gene.2digit.csv",row.names = F,quote = F)

df <- read.csv("20200828/HLAtyping.alle.gene.4digit.csv")
df <- merge(ref[,c(1,2)],df,by.x = "YSample",by.y = "X")
colnames(df)
colnames(df) <- c("YSample","KID","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2","NGS_DRB3.1","NGS_DRB3.2","NGS_DQA1.1", 
                  "NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
write.csv(df,"20200828/HLAtyping.alle.gene.4digit.csv",row.names = F,quote = F)

df <- read.csv("20200828/HLAtyping.alle.gene.6digit.csv")
df <- merge(ref[,c(1,2)],df,by.x = "YSample",by.y = "X")
colnames(df)
colnames(df) <- c("YSample","KID","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2","NGS_DRB3.1","NGS_DRB3.2","NGS_DQA1.1", 
                  "NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
write.csv(df,"20200828/HLAtyping.alle.gene.6digit.csv",row.names = F,quote = F)


##2020
setwd("c:/Users/user/Desktop/KCDC/HLAimputation/HLAtyping/2020/")
ref <- read.csv("HLAtyping.265pairtable.csv",header = T)
head(ref)
df <- read.csv("HLAtyping.alle.gene.2digit.csv")
df <- merge(ref[,c("KBA_ID2","HLAID2")],df,by.x = "HLAID2",by.y = "X")
colnames(df)
colnames(df) <- c("YSample","KID","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2","NGS_DRB3.1","NGS_DRB3.2","NGS_DQA1.1", 
                  "NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
write.csv(df,"HLAtyping.alle.gene.2digit.csv",row.names = F,quote = F)


df <- read.csv("HLAtyping.alle.gene.4digit.csv")
df <- merge(ref[,c("KBA_ID2","HLAID2")],df,by.x = "HLAID2",by.y = "X")
colnames(df)
colnames(df) <- c("YSample","KID","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2","NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2","NGS_DRB3.1","NGS_DRB3.2","NGS_DQA1.1", 
                  "NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2","NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
write.csv(df,"HLAtyping.alle.gene.4digit.csv",row.names = F,quote = F)

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
