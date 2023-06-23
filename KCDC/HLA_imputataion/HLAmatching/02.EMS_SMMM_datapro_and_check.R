##################
## SG
library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)
library(ggrepel)


setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

df <-read.table("Association/c1c2_eplet_weigthSum_association.Result.txt",header = T)
head(df)
c1 <-read_xlsx("Association/c1_forWAS.xlsx")
c2 <-read_xlsx("Association/c2_forWAS.xlsx")


pheno <- read_table("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/00.pheno/Rejection_phenotype_coded_KID_20230323.txt")
ref <- read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt",header = T)




c2_DRDQ <- read_xlsx("Association/c2_eplet_AssoVSS_SMMM.DRDQ.xlsx")
head(c2_DRDQ)
c2_smm_DRDQ <- read_xlsx("Association/KR.KD.HLAimp.foreplet.v2_bCODE_singleMM_scoreCount.xlsx")


colnames(c2)
table(colnames(c2) %in% colnames(c2_smm_DRDQ))
colnames(c2_smm_DRDQ)
c2 %>% select(RecInfo,DonInfo,All_ABV_ClassII,Total_eps,starts_with("Nr")) %>% head()
c2 %>% select(RecInfo,DonInfo,All_ABV_ClassII,Total_eps,starts_with("Nr")) -> c2_nr
### check ALL_ABC_ClassII, Total_eps 
c2 %>% select(RecInfo,DonInfo,All_ABV_ClassII,Total_eps,starts_with("Nr")) %>% 
  #mutate("all_sum" = rowSums(across(starts_with("Nr_all")))) %>% count(all_sum == Total_eps)
  mutate("ABV_sum" = rowSums(across(starts_with("Nr_Ab")))) %>% count(ABV_sum == All_ABV_ClassII)
  
### merged with SMMM (DRB, DQB)
c2 %>% select(RecInfo,DonInfo,All_ABV_ClassII,Total_eps,starts_with("Nr")) %>% 
  select(-colnames(c2_smm_DRDQ)[3:8]) %>% left_join(c2_smm_DRDQ) %>%
  mutate(Total_eps = rowSums(across(starts_with("Nr_all")))) %>%
  mutate(All_ABV_ClassII = rowSums(across(starts_with("Nr_Ab")))) -> c2_smm_nr

head(c2_nr)
head(c2_smm_nr)

c2_nr$theme <- "EMS"
c2_smm_nr$theme <- "SMMS"

c2_nr %>% rbind(c2_smm_nr) %>% select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII) %>%
  pivot_longer(cols = 4:5) %>% #head()
  #mutate(name = factor(name)) %>%
  ggplot(aes(x=value,fill=theme))+
  geom_histogram() + 
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  xlab("# of Eplet count") +
#  theme_ipsum() +
  facet_wrap(~name,scales = "free_x",as.table = FALSE)


