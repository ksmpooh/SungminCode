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

head(pheno)
head(ref)
pheno <- merge(pheno,ref,by.x = "KCHIP_ID",by.y = "KBA_ID",all.x = T)

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

head(c1)
c1 %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>% head()

head(c2_nr)
head(c2_smm_nr)

c2_nr$theme <- "EMS"
c2_smm_nr$theme <- "SMMS"

c2_nr %>% rbind(c2_smm_nr) %>% select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII) %>% #head
  filter(RecInfo %in% pheno$bCODE) %>% 
  pivot_longer(cols = 4:5) %>% 
  group_by(theme,name) %>%
  reframe(q = quantile(value, c(0.25, 0.5, 0.75))) -> c2_quantile

head(c2_quantile)

c2_nr %>% rbind(c2_smm_nr) %>% #head()
  select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII,Nr_allDRB,Nr_AbDRB,Nr_allDQB,Nr_AbDQB) %>% #head
  filter(RecInfo %in% pheno$bCODE) %>% #head()
  pivot_longer(cols = 4:9) %>% #head()
  ggplot(aes(x=value,fill=theme))+
  geom_histogram(binwidth=1) +
  #facet_grid(theme~name,scales = "free_x",as.table = FALSE)  + 
  facet_grid(theme~factor(name,level = c("Total_eps","All_ABV_ClassII","Nr_allDRB","Nr_AbDRB","Nr_allDQB","Nr_AbDQB")),scales = 'free_x') +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  xlab("# of Eplet count") + 
  scale_x_continuous(breaks = pretty_breaks())

c2_nr %>% rbind(c2_smm_nr) %>% select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII) %>% #head
  filter(RecInfo %in% pheno$bCODE) %>% 
  pivot_longer(cols = 4:5) %>% 
  group_by(theme,name) %>%
  mutate(q=ntile(value,4)) %>% 
  group_by(theme,name,q) %>% #count(q)
  mutate(q_m = max(value)) %>% #head()
  mutate(q_m = ifelse(q %in% c(2,4),min(value),q_m))-> a
  #mutate(q_m = ifelse(q %in% c(2,4),))
  count(q)
#-> a

ggplot(a,aes(x=value,fill=theme))+
  geom_histogram(binwidth=1) +
  geom_vline(data = a,aes(xintercept = q_m)) + 
  facet_grid(theme~name,scales = "free_x",as.table = FALSE)  + 
  #facet_grid(name~factor(theme,level = c("Total_eps","All_ABV_ClassII")),scales = "free_x",as.table = FALSE)  + 
  theme(legend.position = 'none',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 13,face = "bold")) + 
  #theme(legend.position = "none",
        #legend.title = element_blank()) +
  xlab("# of Eplet count")
  #  theme_ipsum() +
  #facet_wrap(~name,scales = "free_x",as.table = FALSE)

#### 
c2_nr %>% rbind(c2_smm_nr) %>% select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII) %>%
  filter(RecInfo %in% pheno$bCODE) -> c2_smm_asso

c2_smm_asso %>% hea
c2_smm_asso <- merge(c2_smm_asso,pheno,by.x = "RecInfo",by.y = "bCODE")

result <- NULL
for(i in 4:5){
  temp <- glm(paste("rej_tot ~ ",colnames(c2_smm_asso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_smm_asso %>% filter(theme == "EMS"), family="binomial")
  result <- rbind(result, c("EMS",colnames(c2_smm_asso)[i],as.vector(summary(temp)$coefficients[2,])))
}

for(i in 4:5){
  temp <- glm(paste("rej_tot ~ ",colnames(c2_smm_asso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_smm_asso %>% filter(theme != "EMS"), family="binomial")
  result <- rbind(result, c("SMMS",colnames(c2_smm_asso)[i],as.vector(summary(temp)$coefficients[2,])))
}
result
colnames(result) <- c("theme","ID","Estimate","SE","Z","P")

result <- as.data.frame(result)
result
