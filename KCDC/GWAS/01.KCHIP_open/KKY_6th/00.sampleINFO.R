#IID	FID	CEL	SET	GENDER	BLINDID	COHORT	YEAR	DATATYPE	INFO	BATCH
library(tidyverse)
library(readxl)
library(stringr)
setwd("~/Desktop/KCDC/KKY/00.sampleInfo/")

df <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx")

ori <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx",sheet = 2)
Acc <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx",sheet = 3)
Rep <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과_bysm.xlsx",sheet = 4)

head(ori)
ori %>% select(id,probeset_id,set,gender) %>% 
  mutate(FID = id,BLINDID = "-",COHORT = "KNHANES",YEAR = "2020",DATATYPE="KORV1.1",INFO = "-",BATCH ="KNHANES_6th") %>%
  rename('IID' = id,'CEL' = probeset_id,"SET" = set,"GENDER" = gender) %>% 
  select(IID,FID,CEL,SET,GENDER,BLINDID,COHORT,YEAR,DATATYPE,INFO,BATCH) %>%
  head()


df <- read_excel("2020년도 한국인칩을 이용한 국민건강영양조사 참여자 유전체 정보 생산6136건 결과.xlsx",sheet = 3)
head(df);dim(df)


df <- read.t

