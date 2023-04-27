library(tidyverse)
library(readxl)

setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/NGS(HLAtyping???????????????)/")
setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/NGS(HLAtyping???????????????)/")
df <-read_excel("HLA.typing.Final.result_modify_20211216.xlsx")
head(df)
dim(df)


#FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1
#/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping/all/HLA.type.result.8genes.merged(2019.2020).csv

#df <- read.csv("/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping/all/HLA.type.result.8genes.merged(2019.2020).csv")

df <- read.csv("/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/NGS(HLAtyping???????????????)/HLA.type.result.11genes.merged(2019.2020)_529sample.csv")
colnames(df)
head(df)
head(df[1:5,c("NGS_DPA1.1","NGS_DPA1.2")])
df %>% select(-NGS_DRB3.1,-NGS_DRB3.2,-NGS_DRB4.1,-NGS_DRB4.2,-NGS_DRB5.1,-NGS_DRB5.2) -> df
tcol <- grep("NGS_",colnames(df))

df[,tcol] <- lapply(df[tcol], function(x)
  str_split_fixed(x,"/",2)[,1]
  )
head(df)
df %>% filter(is.na(NGS_A.2))
df %>% filter(NGS_A.2 =="")
df[ df == "." ] <- NA
df[ df == "" ] <- NA
#df[is.na(df) ] <- NA



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


#writexl::write_xlsx(df,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/HLA.typing.Final.result_529sample.xlsx")


#ngs.6d[(!is.na(ngs.6d$DRB1.1)) & (is.na(ngs.6d$DRB1.2)),]$DRB1.2 
#<- ngs.6d[(!is.na(ngs.6d$DRB1.1)) & (is.na(ngs.6d$DRB1.2)),]$DRB1.1


#FID,IID,pID,mID,SEX,PHENO,A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1

df$FID <- df$ID
df$IID <- df$ID
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
write.table(df,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLA.type.result.8genes.merged_529sample_forMAKEreference.txt",col.names = T,row.names = F,quote = F,sep="\t")


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

write.table(df,"/Users/ksmpooh/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLA.type.result.8genes.merged.4digit_529sample_forMAKEreference.txt",col.names = T,row.names = F,quote = F,sep="\t")




############################# 20230306 HLA typing data  merge
######2019+ 2020
library(tidyverse)
library(readxl)
library(stringr)

setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/NGS(HLAtyping?????§á????????)/")
df1 <- read_excel("HLA.typing.Final.result_modify_20211216.xlsx",sheet = 1)
df2 <- read_excel("HLA.typing.Final.result_modify_20211216.xlsx",sheet = 2)
ref <- readxl::read_xls("~/Desktop/KCDC/?????µá???????¡á?«á???????¼á????§á??/?????¡á?¼á????µá????µá????µá??/HLAtyping.265pairtable_with(QC.typing)_20211028.xls")
head(df1)
head(df2)
head(ref)
ref$HLAID.2019 <- str_replace_all(ref$HLAID.2019,"H","CDC")
head(ref)

ref %>% select(HLAID.2019,KBA_ID.2019) %>% rename("Sample" = HLAID.2019,"ID"=KBA_ID.2019) %>%
  inner_join(df1) %>% filter(Sample != 'CDC015') -> out1

ref %>% select(HLAID.2020,KBA_ID.2020) %>% rename("Sample" = HLAID.2020,"ID"=KBA_ID.2020) %>%
  inner_join(df2) -> out2
colnames(out1)
colnames(out2)
out <- rbind(out1,out2)
#write.csv(out,"HLA.type.result.8genes.merged(2019.2020)_529sample.csv")

#### indexing 10 fold 
#setwd("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/2023/")
ref1 <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/HLA.type.result.8genes.merged.4digit_forMAKEreference.txt",header = T) #%>% select(FID)
head(ref)
head(ref1)

c2019 <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/KBA_QC/2019.list.txt")
c2020 <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/KBA_QC/2020.list.txt")

c2019
c2020

head(ref)
ref %>% filter(KBA_ID.2019 %in% c2019$V1) %>% filter(KBA_ID.2020 %in% c2020$V1)

ref %>% filter(KBA_ID.2019 %in% c2020$V1)

qcin <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/KBA_QC/QCinfor.HLAref.txt")
qcout_2020 <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/KBA_QC/QCoutlist.2020.txt")
head(qcin)
head(qcout)
# pair OK
a <- rep(c(1:10), times=100)
ref %>% filter( HLAID.2019 != 'H015') %>%
  count(typing2019,typing2020,QC.2020) #%>%

# A tibble: 3 x 4
'
typing2019 typing2020 QC.2020     n
<dbl>      <dbl> <chr>   <int>
1          0          2 Pass       10
2          1          0 Fail        3
3          1          1 Pass      251  
'
# select : 
#typing 2020 : 259 + KBA quality good 1), 
#typing 2019 : 254 + KBA quality good 6)
head(ref$HLAID.2019)
ref %>% mutate(typing2019 = ifelse(KBA_ID.2019 %in% qcin$V1,1,typing2019)) %>%
  mutate(typing2020 = ifelse(KBA_ID.2020 %in% qcin$V1,1,typing2020)) %>%
  mutate(typing2020 = ifelse(KBA_ID.2020 %in% qcout_2020$V1,3,typing2020)) %>% 
  filter( HLAID.2019 != 'CDC015') %>%
  count(typing2019,typing2020,QC.2020)
'
typing2019 typing2020 QC.2020     n
<dbl>      <dbl> <chr>   <int>
1          0          2 Pass        2 2020
2          0          3 Pass        2 x
3          1          0 Fail        3 pair
4          1          1 Pass      252 pair (+CDC015)
5          1          2 Pass        3 pair
6          1          3 Pass        3 2019
'

ref %>% mutate(typing2019 = ifelse(KBA_ID.2019 %in% qcin$V1,1,typing2019)) %>%
  mutate(typing2020 = ifelse(KBA_ID.2020 %in% qcin$V1,1,typing2020)) %>%
  mutate(typing2020 = ifelse(KBA_ID.2020 %in% qcout_2020$V1,3,typing2020)) -> ref

head(ref)
# 258 pair = 516 sample
a <- rep(c(1:5), times=100)
ref %>%
  filter( HLAID.2019 != 'CDC015') %>%
  #count(typing2019,typing2020,QC.2020)
  filter(typing2019 != 0,typing2020 !=3) %>% #dim()
  select(KBA_ID.2019,KBA_ID.2020) %>%
  mutate('group' = a[1:257]) -> p1


head(p1)
p2 <- p1 %>% select(KBA_ID.2019,group) %>% rename("KBAID" = KBA_ID.2019)
p1 <- p1 %>% select(KBA_ID.2020,group) %>% rename("KBAID" = KBA_ID.2020)
head(p1)
head(p2)

p1 %>% full_join(p2) %>% count(group)


# 4 sample 
a <- rep(c(3:5), times=5)
ref %>% 
  filter((typing2020 == 2 & typing2019 == 0) | HLAID.2019 == 'CDC015') %>% #count(typing2019,typing2020)
  select(KBA_ID.2020) %>% rename("KBAID" = KBA_ID.2020) %>%
  mutate('group' = a[1:3]) -> p3

a <- rep(c(3:5), times=5)
ref %>% 
  filter((typing2020 == 3 & typing2019 == 1)) %>% #count(typing2019,typing2020)
  select(KBA_ID.2019) %>% rename("KBAID" = KBA_ID.2019) %>%
  mutate('group' = a[1:3]) -> p4

head(p3)
head(p4)


table(a$group)
a <- p1 %>% full_join(p2) %>% full_join(p3) %>% full_join(p4) #%>% count(group)#%>% count(group)
head(a)
write.table(a,"~/Desktop/KCDC/HLAimputation/MakeReferencePanel/Final_520sample.index.txt",col.names = T,row.names = F,quote = F,sep = "\t")


1.1 == 1 + 0.1
