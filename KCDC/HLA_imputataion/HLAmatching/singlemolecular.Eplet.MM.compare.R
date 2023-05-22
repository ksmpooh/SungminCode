#### HLA eplet single molecular

library(tidyverse)
library(readxl)
library(writexl)
library(readxlsb)


setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")


c2_1 <- read_xlsx("20230508_epletv2/KR.KD.HLAimp.foreplet.v2_bCODE.xlsx",sheet = 3,col_names = F)
c2_2 <- read_xlsx("20230508_epletv2/KR.KD.HLAimp.foreplet.v2_bCODE.xlsx",sheet = 4,col_names = F)

head(c2_1)
head(c2_2)

c2 <- rbind(c2_1,c2_2)
head(c2)
colnames(c2) <- c("KR","Rec1stDRB","Rec2ndDRB","Rec1stDRW","Rec2ndDRW","Rec1stDQB","Rec2ndDQB","Rec1stDQA",
                    "Rec2ndDQA","Rec1stDPB","Rec2ndDPB","Rec1stDPA","Rec2ndDPA","KD","D1stDRB","D2ndDRB","D1stDRW",
                    "D2ndDRW","D1stDQB","D2ndDQB","D1stDQA","D2ndDQA","D1stDPB","D2ndDPB","D1stDPA","D2ndDPA")

#KR_bCODE,HLA_DRB1.1.x,HLA_DRB1.2.x,x,x1,HLA_DQB1.1.x,HLA_DQB1.2.x,HLA_DQA1.1.x,HLA_DQA1.2.x,HLA_DPB1.1.x,HLA_DPB1.2.x,HLA_DPA1.1.x,HLA_DPA1.2.x,KD_bCODE,HLA_DRB1.1.y,HLA_DRB1.2.y,x2,x3,HLA_DQB1.1.y,HLA_DQB1.2.y,HLA_DQA1.1.y,HLA_DQA1.2.y,HLA_DPB1.1.y,HLA_DPB1.2.y,HLA_DPA1.1.y,HLA_DPA1.2.y
colnames(c2)

c2_h1 <- c2[,c(1,2,3,4,4,6,7,4,4,4,4,4,4,14,15,15,4,4,19,19)]
c2_h2 <- c2[,c(1,2,3,4,4,6,7,4,4,4,4,4,4,14,16,16,4,4,20,20)]

head(c2_h1)

#dataset_names <- list("classII_g1_h1" = c2_h1[1:1000,], "classII_g2_h1" = c2_h1[1001:nrow(c2_h1),],"classII_g1_h2" = c2_h2[1:1000,], "classII_g2_h2" = c2_h2[1001:nrow(c2_h2),])
#writexl::write_xlsx(dataset_names,"20230508_epletv2/single_molecular_MM/KR.KD.HLAimp.foreplet.v2_bCODE_singleMM.xlsx",col_names = F)


##### eplet score comparison 
#### after input HLA match maker

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230508_epletv2/")

c2_g1_h1 <- read_xlsb("single_molecular_MM/DRDQDP_Eplet_Matching_3.1_g1_h1.xlsb",sheet = 4,skip = 4)
c2_g1_h2 <- read_xlsb("single_molecular_MM/DRDQDP_Eplet_Matching_3.1_g1_h2.xlsb",sheet = 4,skip = 4)
c2_g2_h12 <- read_xlsb("single_molecular_MM/DRDQDP_Eplet_Matching_3.1_g2_h12.xlsb",sheet = 4,skip = 4)


head(c2_g1_h1)
head(c2_g1_h2)
head(c2_g2_h12)

c2_g1_h1 <- c2_g1_h1 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")
c2_g1_h2 <- c2_g1_h2 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")

c2_g2_h12 <- c2_g2_h12 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")

dim(c2_g2_h12)

c2_h1 <- rbind(c2_g1_h1,c2_g2_h12[1:147,])
c2_h2 <- rbind(c2_g1_h2,c2_g2_h12[148:294,])

head(c2_h2)
head(c2_h1)

'
c2_header <- read_xlsb("DRDQDP_Eplet_Matching_3.1.xlsb",sheet = 4,skip = 1)
#c2_header[1,] %>% t()

c2_header[1,] %>% write_xlsx("c2_header_extra.xlsx")
# 이거 후에 수작업
'

c2_hh <- read_excel("c2_header.xlsx",sheet = 4)
c2_hh
c2_hh$V5

#c2_hh %>% t() %>% as.data.frame() -> c2_header_index 

colnames(c2_h1) <- seq(1:ncol(c2_h1))
colnames(c2_h2) <- seq(1:ncol(c2_h2))

c2_hh$V5[grepl("DR",c2_hh$V5)]
grepl("DR",c2_hh$V5)

head(c2_h1)

c2_h1 %>% select(c2_hh$V4) -> c2_h1
c2_h2 %>% select(c2_hh$V4) -> c2_h2

colnames(c2_h1) <- c2_hh$V5
colnames(c2_h2) <- c2_hh$V5


head(c2_h1)
head(c2_h2)

c2_hh$V5

c2_h1_sm <- c2_h1[,c("RecInfo","DonInfo","Nr_allDRB","Nr_AbDRB","Nr_otDRB", "Nr_allDQB","Nr_AbDQB","Nr_otDQB")]
c2_h2_sm <- c2_h1[,c("RecInfo","DonInfo","Nr_allDRB","Nr_AbDRB","Nr_otDRB", "Nr_allDQB","Nr_AbDQB","Nr_otDQB")]
head(c2_h1_sm)
head(c2_h2_sm)

c2_h1_sm %>% #merge(c2_h2_sm,by = "RecInfo") %>% #count(DonInfo.x == DonInfo.y)
  merge(c2_h2_sm,by = c("RecInfo","DonInfo")) %>% #head()
  mutate("Nr_allDRB" = ifelse(Nr_allDRB.x > Nr_allDRB.y,Nr_allDRB.x,Nr_allDRB.y)) %>% #head()
  mutate("Nr_AbDRB" = ifelse(Nr_AbDRB.x > Nr_AbDRB.y,Nr_AbDRB.x, Nr_AbDRB.y)) %>% #head()
  mutate("Nr_otDRB" = ifelse(Nr_otDRB.x > Nr_otDRB.y,Nr_otDRB.x,Nr_otDRB.y)) %>% #head()
  mutate("Nr_allDQB" = ifelse(Nr_allDQB.x > Nr_allDQB.y,Nr_allDQB.x,Nr_allDQB.y)) %>% #head()
  mutate("Nr_AbDQB" = ifelse(Nr_AbDQB.x > Nr_AbDQB.y,Nr_AbDQB.x,Nr_AbDQB.y)) %>% #head()
  mutate("Nr_otDQB" = ifelse(Nr_otDQB.x > Nr_otDQB.y,Nr_otDQB.x,Nr_otDQB.y)) %>% #head()
  select(RecInfo,DonInfo,Nr_allDRB,Nr_AbDRB,Nr_otDRB,Nr_allDQB,Nr_AbDQB,Nr_otDQB) -> c2_mm_forASO


write_xlsx(c2_mm_forASO,"single_molecular_MM/KR.KD.HLAimp.foreplet.v2_bCODE_singleMM_scoreCount.xlsx")

### asso 비교

pheno <- read_table("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/00.pheno/Rejection_phenotype_coded_KID_20230323.txt")
ref <- read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt",header = T)

#c2 <- read_xlsx("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/Association/c2_forWAS.xlsx")
#c2 <- read_table("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/Association/c2_eplet_asso.txt")
c2 <- read_table("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/Association/c2_eplet_asso_withouthHLAms.txt")
c2_mm_forASO <- read_xlsx("single_molecular_MM/KR.KD.HLAimp.foreplet.v2_bCODE_singleMM_scoreCount.xlsx")
head(c2_mm_forASO)
pheno %>% merge(ref,by.x = "KCHIP_ID",by.y="KBA_ID") %>% #head()
  merge(c2_mm_forASO,by.x="bCODE",by.y="RecInfo") -> c2_mm_data


head(c2_mm_data)
table(c2_mm_data$rej_tot)
colnames(c2_data)[43:ncol(c2_data)]
colnames(c2_data)[60:ncol(c2_data)]

colnames(c2_mm_data) %>% as_data_frame() %>% count(value) %>% filter(n==2)
colnames(c2_mm_data)

#colnames(data2)[36:ncol(data2)] <- paste0("AA_",colnames(data2)[36:ncol(data2)])

str(c2_mm_data)
'
for (i in 18:ncol(c2_mm_data)) {
  c2_data[,i] <- as.numeric(c2_mm_data[,i])
}
'


result <- NULL
for(i in 18:ncol(c2_mm_data)){
  #table(c2_data[,i])
  temp <- glm(paste("rej_tot ~ ",colnames(c2_mm_data)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_mm_data, family="binomial")
  result <- rbind(result, c(colnames(c2_mm_data)[i],as.vector(summary(temp)$coefficients[2,])))
}

head(result)

colnames(result) <- c("ID","Estimate","SE","Z","P")
c2_mm_result <- result %>% as.data.frame()

head(c2_mm_result)
c2_mm_result$theme <- "SMMM"
c2$theme <- "EMS"

head(c2)
head(c2_mm_result)

c2_mm_result %>% rbind(c2[c2$ID %in% c2_mm_result$ID,]) %>% #head()
  pivot_longer(cols = 2:5,names_to = "Stat.",values_to = "Value") %>% #head()
  pivot_wider(names_from = "theme",values_from = "Value") -> tt

#pivot_wider(names_from = "theme",values_from = "Value")c2_mm_result %>% rbind(c2[c2$ID %in% c2_mm_result$ID,]) %>% head()
  
c2_mm_result %>% rbind(c2[c2$ID %in% c2_mm_result$ID,]) %>% #head()
  arrange(ID,theme) %>% writexl::write_xlsx("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/Association/c2_eplet_AssoVSS_SMMM.DRDQ.xlsx")
    
head(tt)


tt %>% ggplot(aes(x=SMMM,y=EMS,color=Stat.)) +
  geom_point() +
  facet_grid(~ID)

head(pheno)
dim(pheno)
dim(ref)
head(c2_mm_data)
dim(c2_mm_data) 
table(c2_mm_data$rej_tot)
