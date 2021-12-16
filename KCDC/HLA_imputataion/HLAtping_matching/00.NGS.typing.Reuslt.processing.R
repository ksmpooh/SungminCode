##### HLA imputation 결과 정리

library(readxl)
library(writexl)
library(tidyverse)
library(stringr)


##ref
ref <- read_excel("~/Desktop/KCDC/일반용역/HLAtyping.265pairtable_with(QC.typing)_20211028.xls")%>%
  select(KBA_ID.2019,oriID.2019,HLAID.2019,KBA_ID.2020,oriID.2020,HLAID.2020) %>% as.data.frame() 
colnames(ref)                  
head(ref)
ref$HLAID.2019 <- str_replace_all(ref$HLAID.2019,"H","CDC")
ref <- ref %>% filter(HLAID.2019 != "CDC015")


df <- read_excel("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/HLA.typing.Final.result_modify_20211216.xlsx",sheet = 1)
head(df)
df1 <- read_excel("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/HLA.typing.Final.result_modify_20211216.xlsx",sheet = 2)
head(df1)

df <- rbind(df,df1)
colnames(df)
head(df)

ngs <-df
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


ngs.2d <- subset(ngs,select =c("Sample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DRB3.1","DRB3.2","DQA1.1",
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



ngs.4d <- subset(ngs,select =c("Sample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DRB3.1","DRB3.2","DQA1.1",
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

ngs.6d <- subset(ngs,select =c("Sample","A.1","A.2","B.1","B.2","C.1","C.2","DRB1.1","DRB1.2","DRB3.1","DRB3.2","DQA1.1",
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


head(ref)

a <- ref %>% select(KBA_ID.2019,oriID.2019,HLAID.2019) %>% rename(KBA_ID = "KBA_ID.2019",OriID = 'oriID.2019',HLAID = 'HLAID.2019')
b <- ref %>% select(KBA_ID.2020,oriID.2020,HLAID.2020) %>% rename(KBA_ID = "KBA_ID.2020",OriID = 'oriID.2020',HLAID = 'HLAID.2020')

ref1 <- rbind(a,b)
head(ref1)
head(ngs.2d)

ngs.2d <- merge(ref1,ngs.2d,by.x = "HLAID",by.y="Sample")
ngs.4d <- merge(ref1,ngs.4d,by.x = "HLAID",by.y="Sample")
ngs.6d <- merge(ref1,ngs.6d,by.x = "HLAID",by.y="Sample")

write.table(ngs.2d,"~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/HLA.typing_2digit.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(ngs.4d,"~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/HLA.typing_4digit.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(ngs.6d,"~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/HLA.typing_6digit.txt",col.names = T,row.names = F,quote = F,sep = "\t")


