## HLA typing result processing

setwd("c:/Users/user/Desktop/KCDC/HLAimputation/HLAtyping/all/")
ori <- read.csv("HLA.type.result.8genes.merged(2019.2020).csv")
head(ori)
colnames(ori) <- c("KID","YSample","A_1","A_2","B_1","B_2","C_1","C_2","DRB1_1","DRB1_2"
                   ,"DQA1_1","DQA1_2","DQB1_1","DQB1_2"
                   ,"DPA1_1","DPA1_2","DPB1_1","DPB1_2")
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

ngs$DQA1.1 <- str_split_fixed(str_split_fixed(ngs$DQA1_1,":",2)[,1],"\\*",2)[,2]
ngs$DQA1.2 <- str_split_fixed(str_split_fixed(ngs$DQA1_2,":",2)[,1],"\\*",2)[,2]
ngs$DQB1.1 <- str_split_fixed(str_split_fixed(ngs$DQB1_1,":",2)[,1],"\\*",2)[,2]
ngs$DQB1.2 <- str_split_fixed(str_split_fixed(ngs$DQB1_2,":",2)[,1],"\\*",2)[,2]

ngs$DPA1.1 <- str_split_fixed(str_split_fixed(ngs$DPA1_1,":",2)[,1],"\\*",2)[,2]
ngs$DPA1.2 <- str_split_fixed(str_split_fixed(ngs$DPA1_2,":",2)[,1],"\\*",2)[,2]
ngs$DPB1.1 <- str_split_fixed(str_split_fixed(ngs$DPB1_1,":",2)[,1],"\\*",2)[,2]
ngs$DPB1.2 <- str_split_fixed(str_split_fixed(ngs$DPB1_2,":",2)[,1],"\\*",2)[,2]


colnames(ngs)


ngs.2d <- subset(ngs,select =c("KID","YSample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.2d)
features <- colnames(ngs.2d)[3:ncol(ngs.2d)]
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

head(ngs.2d)

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

ngs$DQA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_1,":",3)[,1],str_split_fixed(ngs$DQA1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DQA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_2,":",3)[,1],str_split_fixed(ngs$DQA1_2,":",3)[,2]),"\\*",2)[,2]
ngs$DQB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_1,":",3)[,1],str_split_fixed(ngs$DQB1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DQB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_2,":",3)[,1],str_split_fixed(ngs$DQB1_2,":",3)[,2]),"\\*",2)[,2]

ngs$DPA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_1,":",3)[,1],str_split_fixed(ngs$DPA1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DPA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_2,":",3)[,1],str_split_fixed(ngs$DPA1_2,":",3)[,2]),"\\*",2)[,2]
ngs$DPB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_1,":",3)[,1],str_split_fixed(ngs$DPB1_1,":",3)[,2]),"\\*",2)[,2]
ngs$DPB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_2,":",3)[,1],str_split_fixed(ngs$DPB1_2,":",3)[,2]),"\\*",2)[,2]



ngs.4d <- subset(ngs,select =c("KID","YSample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.4d)
features <- colnames(ngs.4d)[3:ncol(ngs.4d)]
features
for (i in features){
  #ngs.4d[,i] <- strtoi(ngs.4d[,i])
  ngs.4d[,i] <- as.integer(ngs.4d[,i])
}

ngs.4d[(!is.na(ngs.4d$DRB1.1)) & (is.na(ngs.4d$DRB1.2)),]$DRB1.2 <- ngs.4d[(!is.na(ngs.4d$DRB1.1)) & (is.na(ngs.4d$DRB1.2)),]$DRB1.1
ngs.4d[(is.na(ngs.4d$DRB1.1)) & (!is.na(ngs.4d$DRB1.2)),]$DRB1.1 <- ngs.4d[(is.na(ngs.4d$DRB1.1)) & (!is.na(ngs.4d$DRB1.2)),]$DRB1.2

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

ngs$DQB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_1,":",4)[,1],str_split_fixed(ngs$DQB1_1,":",4)[,2],str_split_fixed(ngs$DQB1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DQB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQB1_2,":",4)[,1],str_split_fixed(ngs$DQB1_2,":",4)[,2],str_split_fixed(ngs$DQB1_2,":",4)[,3]),"\\*",2)[,2]
ngs$DQA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_1,":",4)[,1],str_split_fixed(ngs$DQA1_1,":",4)[,2],str_split_fixed(ngs$DQA1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DQA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DQA1_2,":",4)[,1],str_split_fixed(ngs$DQA1_2,":",4)[,2],str_split_fixed(ngs$DQA1_2,":",4)[,3]),"\\*",2)[,2]

ngs$DPA1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_1,":",4)[,1],str_split_fixed(ngs$DPA1_1,":",4)[,2],str_split_fixed(ngs$DPA1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DPA1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPA1_2,":",4)[,1],str_split_fixed(ngs$DPA1_2,":",4)[,2],str_split_fixed(ngs$DPA1_2,":",4)[,3]),"\\*",2)[,2]
ngs$DPB1.1 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_1,":",4)[,1],str_split_fixed(ngs$DPB1_1,":",4)[,2],str_split_fixed(ngs$DPB1_1,":",4)[,3]),"\\*",2)[,2]
ngs$DPB1.2 <- str_split_fixed(paste0(str_split_fixed(ngs$DPB1_2,":",4)[,1],str_split_fixed(ngs$DPB1_2,":",4)[,2],str_split_fixed(ngs$DPB1_2,":",4)[,3]),"\\*",2)[,2]

ngs.6d <- subset(ngs,select =c("KID","YSample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DQA1.1",
                               "DQA1.2","DQB1.1","DQB1.2","DPA1.1","DPA1.2","DPB1.1","DPB1.2"))
str(ngs.6d)
features <- colnames(ngs.6d)[3:ncol(ngs.6d)]
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

colnames(ngs.2d) <- c("KID","YSample","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2"
                      ,"NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2"
                      ,"NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2"
                    ,"NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
colnames(ngs.4d) <- c("KID","YSample","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2"
                      ,"NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2"
                      ,"NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2"
                      ,"NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")
colnames(ngs.6d) <- c("KID","YSample","NGS_A.1","NGS_A.2","NGS_B.1","NGS_B.2"
                      ,"NGS_C.1","NGS_C.2","NGS_DRB1.1","NGS_DRB1.2"
                      ,"NGS_DQA1.1","NGS_DQA1.2","NGS_DQB1.1","NGS_DQB1.2"
                      ,"NGS_DPA1.1","NGS_DPA1.2","NGS_DPB1.1","NGS_DPB1.2")


#2019 + 2020 
write.csv(ngs.2d,"HLAtyping.alle.gene.2digit_2019.with.2020.csv",row.names = F,quote = F)
write.csv(ngs.4d,"HLAtyping.alle.gene.4digit_2019.with.2020.csv",row.names = F,quote = F)
write.csv(ngs.6d,"HLAtyping.alle.gene.6digit_2019.with.2020.csv",row.names = F,quote = F)


