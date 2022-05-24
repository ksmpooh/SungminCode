#IID	FID	CEL	SET	GENDER	BLINDID	COHORT	YEAR	DATATYPE	INFO	BATCH
library(tidyverse)
library(readxl)
library(stringr)
setwd("~/Desktop/KCDC/KKY/00.sampleInfo/")

df <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx")

ori <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx",sheet = 2)
Rep <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx",sheet = 3)
Acc <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx",sheet = 4)

ref <- read.table("2020.all.teragen.KKY.info_20220311.txt",header = T)
ref$BATCH <- paste0("KNHANES_",ref$th)
head(ref)

mzlist1 <- read.table("2020_tera/KBA_Blind_matchingTable.txt")
mzlist2 <- read.table("2020_tera/KBA_TPMRA_matchingtable_v2.txt")

head(mzlist)
str_split_fixed(str_split_fixed(mzlist$V1,"_",6)[,6],"\\.",2)[,1]

mzlist1 <- mzlist1 %>% select(V1,V3) %>% 
  mutate(IID=str_split_fixed(str_split_fixed(V1,"_",6)[,6],"\\.",2)[,1],BLINDID=str_split_fixed(str_split_fixed(V3,"_",6)[,6],"\\.",2)[,1]) %>%
  rename("KCHIP_CEL" = V1, "BLIND_CEL" = V3)

mzlist2 <- mzlist2 %>% select(V1,V3) %>% 
  mutate(IID=str_split_fixed(str_split_fixed(V1,"_",6)[,6],"\\.",2)[,1],BLINDID=str_split_fixed(str_split_fixed(V3,"_",5)[,5],"\\.",2)[,1]) %>%
  rename("KCHIP_CEL" = V1, "BLIND_CEL" = V3)

head(mzlist1)
head(mzlist2)
mzlist <- rbind(mzlist1,mzlist2)

head(ori)
a <- ori %>% select(id,probeset_id,set,gender) %>% 
  mutate(FID = id,COHORT = "KNHANES",YEAR = "2020",DATATYPE="KORV1.1",INFO = "-") %>%
  rename('IID' = id,'CEL' = probeset_id,"SET" = set,"GENDER" = gender) %>% 
  merge(ref[,c(1,3)],by.x = "FID",by.y = "id") %>%
  full_join(mzlist[,c("IID","BLINDID")]) %>%
  mutate(BLINDID = ifelse(is.na(BLINDID),"-",BLINDID)) %>%
  select(IID,FID,CEL,SET,GENDER,BLINDID,COHORT,YEAR,DATATYPE,INFO,BATCH) %>%
  mutate(SET = paste0("DL",str_replace(SET,"'","")))

table(a$BLINDID)
head(a)
a %>% filter(FID == "NIH20C2841404")
mzlist %>% filter(IID == "NIH20C2841404")
head(mzlist)

head(mzlist)
colnames(mzlist)[3:4] <- c("BLINDID","IID")
b <-Rep %>% select(id,probeset_id,set,gender) %>% 
  mutate(FID = id,COHORT = "KNHANES",YEAR = "2020",DATATYPE="KORV1.1",INFO = "Reproducibility") %>%
  rename('IID' = id,'CEL' = probeset_id,"SET" = set,"GENDER" = gender) %>% 
  inner_join(mzlist[,c("IID","BLINDID")]) %>%
  merge(ref[,c(1,3)],by.x = "BLINDID",by.y = "id") %>% 
  select(IID,FID,CEL,SET,GENDER,BLINDID,COHORT,YEAR,DATATYPE,INFO,BATCH) %>%
  mutate(SET = paste0("DL",str_replace(SET,"'","")))

head(b)

c <- Acc %>% select(id,probeset_id,set,apt_probeset_genotype_gender) %>% 
  mutate(FID = id,COHORT = "KNHANES",YEAR = "2020",DATATYPE="TPMRA",INFO = "Accuracy") %>%
  rename('IID' = id,'CEL' = probeset_id,"SET" = set,"GENDER" = apt_probeset_genotype_gender) %>% 
  inner_join(mzlist[,c("IID","BLINDID")]) %>%
  merge(ref[,c(1,3)],by.x = "BLINDID",by.y = "id") %>% 
  select(IID,FID,CEL,SET,GENDER,BLINDID,COHORT,YEAR,DATATYPE,INFO,BATCH) %>%
  mutate(SET = paste0("DL",str_replace(SET,"'","")))


head(c)
head(mzlist)
head(ref)
df <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과.xlsx",sheet = 3)
head(df);dim(df)
head(a)
df <- read.table("2020_tera/gender_Mismatch.txt",header = F)
gender_ref <- read.table("2020_tera/gender.mismatch.id.txt",header =F) %>% 
  mutate(BATCH=ifelse(V1 == "2013" || V1 =="2014" ||V1 == "2015","KNHANES_6th","KNHANES_7th")) %>%
  rename("FID" = V2)
head(gender_ref)


head(df)
d<-df %>% mutate(COHORT = "KNHANES",YEAR = "2020",DATATYPE="KORV1.1",INFO = "gender_discrepancy_fail",BLINDID="-",GENDER="-") %>%
  mutate(FID = str_split_fixed(str_split_fixed(V1,"_",6)[,6],"\\.",2)[,1]) %>%
  mutate(IID = FID) %>%
  rename("CEL" = V1) %>%
  inner_join(gender_ref[,c(2,3)]) %>%
  mutate(SET = str_split_fixed(CEL,"_",6)[,4]) %>%
  mutate(SET = paste0("DL",str_replace(SET,"'",""))) %>%
  select(IID,FID,CEL,SET,GENDER,BLINDID,COHORT,YEAR,DATATYPE,INFO,BATCH)
  
out <- rbind(a,b)
out <- rbind(out,c)    
out <- rbind(out,d)
dim(out)
head(out)


library(writexl)
write_xlsx(out,"2020_tera/2020.PRO.FINAL.info.xlsx")
