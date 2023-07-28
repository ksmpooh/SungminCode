### acute rejection new pheno  20230705
library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)
library(ggrepel)
library(scales)

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

pheno <- read_xlsx("00.pheno/전희중교수님 주신 파일_20221109_KOTRY_신장_1912 검체 매칭 2395_rdm(rejection)추가변수.xlsx")
#ref <-read.table("00.pheno/Rejection_phenotype_coded_KID_20230323.txt",header = T)
pair <- read.table("00.pheno/02.KR.KD.immune.cell.co-signal_targetGene.alleleMatching01.Score_Sum.txt",header = T)
ref <-read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt",header = T)
head(pair)[1:5]
pair %>% merge(ref,by.x = "KBA_ID.KR",by.y="KBA_ID") %>% select(KBA_ID.KR,bCODE) -> pair


head(ref)
dim(ref)
head(pheno)
colnames(pheno)
colnames(pheno)[grep("group",colnames(pheno))]

colnames(pheno)[grep("hla",colnames(pheno))]

pheno %>% count(group...11 == group...10268)
pheno %>% count(group...11 == group...10374)
pheno %>% count(cvd == CVD)
pheno %>% count(cmv_igg_reci...10279 == cmv_igg_reci...232)
pheno %>% count(hbsag_reci...236 == hbsag_reci...10280)
pheno %>% count(CYCDR_B...329 == CYCDR_B...10286)
pheno %>% select(`인체자원:자원bCODE...2`,`인체자원:자원bCODE...4`,rej_tot,date_rej_tot,date_kt,AGE,`SEX...14`,dm,cvd,cmv_igg_reci...10279,hbsag_reci...236,hcv_ab_reci...239,ind_atg...10282,D_AGE,D_SEX...32,dgf...521,TACDR_B...10284,CYCDR_B...329,desen...77,hla_ms_total,`group...11`,hla_a1_reci...248,hla_a2_reci...249,hla_b1_reci...250,hla_b2_reci...251,hla_dr1_reci...252,hla_dr2_reci...253) -> new_pheno
head(new_pheno)
colnames(new_pheno)
colnames(new_pheno) <- c("RecInfo","DonInfo","rej_tot","date_rej_tot","date_kt","AGE","SEX","dm","cvd","cmv_igg_reci","hbsag_reci","hcv_ab_reci",
                         "ind_atg","D_AGE","D_SEX","dgf","TACDR","CYCDR","desen","hla_ms_total","group",
                         "hla_a1_reci","hla_a2_reci","hla_b1_reci","hla_b2_reci","hla_dr1_reci","hla_dr2_reci")

new_pheno %>% filter(RecInfo %in% pair$bCODE) -> a
head(a)
colnames(a)["dm"]3
table(a[,colnames(a)[i]])
dim(a)
'''
Y값: rej_tot
Time값: date_rej_tot - date_kt => month 수치로 변환 필요
Covariates (김형우 교수 메일, 2023/5/29)
Covariate 목록: AGE, SEX, dm, cvd, cmv_igg_reci, hbsag_reci, hcv_ab_reci, ind_atg, D_AGE, D_SEX, dgf, TACDR, CYCDR, desen, hla_ms_tot, group
AGE, D_AGE, hla_tot_ms 제외 모두 binary (0, 1 형태인지 체크 필요), group의 경우는 1,2 형태인지 체크 필요
빈칸으로 되어있는 것은 NA 처리
AGE, D_AGE => '세' 제외하고 숫자로
SEX, D_SEX => Male 0, Female 1로 변경
'''


a[,c("rej_tot","dm","cvd","cmv_igg_reci","hbsag_reci","hcv_ab_reci","ind_atg","dgf","TACDR","CYCDR","desen")] %>% 
  pivot_longer(1:11) %>% count(name,value) -> a1
for (i in c("rej_tot","dm","cvd","cmv_igg_reci","hbsag_reci","hcv_ab_reci","ind_atg","dgf","TACDR","CYCDR","desen")) {
  print(i)
  
}


new_pheno %>% filter(RecInfo %in% pair$bCODE) %>% #head() #count(rej_tot)
  mutate('rej_time_day' = date_rej_tot - date_kt) %>% #select(rej_time) ->a
  mutate(across("AGE",str_replace,"세","")) %>% 
  mutate(across("D_AGE",str_replace,"세","")) %>% #select(D_AGE)
  mutate(AGE = as.numeric(AGE),D_AGE = as.numeric(D_AGE)) %>% #select(SEX)
  mutate(SEX = ifelse(SEX == "Male",0,1)) %>%
  mutate(D_SEX = ifelse(D_SEX == "Male",0,1)) %>% #select(SEX)
  merge(pair,by.x = "RecInfo",by.y="bCODE") %>%
  #writexl::write_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso.xlsx")
  writexl::write_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso.xlsx")
  
df <- read_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso_v2.xlsx")  
df <- read_xlsx("20221109_KOTRY_selectFeature_forAuteRejectionAsso_v3.xlsx")
head(df)  
dim(df)

colnames(df)
df %>% mutate('rej_time_day' = date_rej_tot - date_kt) %>% #select(rej_time_day) %>% #count(rej_time_day)
  writexl::write_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso_v2.xlsx")
#  glm(paste("rej_tot ~ group+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=., family="binomial") -> a
  
## sample check

dim(new_pheno)

sampleInfo <- read_xlsx("~/Desktop/KCDC/transplantation/00.sampleInfo/ALL(2019and2020).sampleID.xlsx")
head(ref)
dim(ref)
dim(new_pheno)
head(new_pheno)
sampleInfo %>% merge(ref) %>% 
  merge(new_pheno,by.x="bCODE",by.y = 'RecInfo',all.y = T) %>% #head( )
  count(type,Prod)

