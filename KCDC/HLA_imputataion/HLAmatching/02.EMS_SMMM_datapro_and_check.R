##################
## SG
library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)
library(ggrepel)
library(scales)



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
head(c2_smm_DRDQ)

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
#  filter(RecInfo %in% pheno$bCODE) %>% 
  pivot_longer(cols = 4:5) %>% 
  group_by(theme,name) %>%
  reframe(q = quantile(value, c(0.25, 0.5, 0.75))) -> c2_quantile

head(c2_quantile)

c2_nr %>% rbind(c2_smm_nr) %>% #head()
  select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII,Nr_allDRB,Nr_AbDRB,Nr_allDQB,Nr_AbDQB) %>% #head
#  filter(RecInfo %in% pheno$bCODE) %>% #head()
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
#  filter(RecInfo %in% pheno$bCODE) %>% 
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
c2_nr %>% rbind(c2_smm_nr) %>% head()
c2_nr %>% rbind(c2_smm_nr) %>% select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII) %>%
#  filter(RecInfo %in% pheno$bCODE) %>% #head()
  mutate("All_Oth_ClassII" = Total_eps-All_ABV_ClassII) -> c2_smm_asso

c2_smm_asso %>% head()
c2_smm_asso <- merge(c2_smm_asso,pheno,by.x = "RecInfo",by.y = "bCODE")

result <- NULL
for(i in 4:6){
  temp <- glm(paste("rej_tot ~ ",colnames(c2_smm_asso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_smm_asso %>% filter(theme == "EMS"), family="binomial")
  result <- rbind(result, c("EMS",colnames(c2_smm_asso)[i],as.vector(summary(temp)$coefficients[2,])))
}

for(i in 4:6){
  temp <- glm(paste("rej_tot ~ ",colnames(c2_smm_asso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_smm_asso %>% filter(theme != "EMS"), family="binomial")
  result <- rbind(result, c("SMMS",colnames(c2_smm_asso)[i],as.vector(summary(temp)$coefficients[2,])))
}

result_scale <- NULL

result
colnames(result) <- c("theme","ID","Estimate","SE","Z","P")

result <- as.data.frame(result)
class2_total <- result



######### Eplet score data scale
#c2_nr %>% rbind(c2_smm_nr) %>% select(RecInfo,DonInfo,theme,Total_eps,All_ABV_ClassII) %>% #head
head(c1)
c2_nr %>% rbind(c2_smm_nr) %>% 
  #filter(RecInfo %in% pheno$bCODE) %>% #head()
  mutate("All_Oth_ClassII" = Total_eps-All_ABV_ClassII) -> c2_smm_asso

c1 %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) -> c1_asso

head(c1_asso)
colnames(c1_asso) <- c("RecInfo","DonInfo","Total_eps_ClassI","Abv_ClassI","Oth_ClassI")

'''
colnames(c2_smm_asso) <- c("RecInfo","DonInfo","Abv_ClassII","Total_eps_ClassII","allDRB","AbDRB","otDRB","allDQB",
                           "AbDQB","otDQB","allDQA","AbDQA","otDQA","allDPB","AbDPB","otDPB",
                           "allDPA","AbDPA","otDPA","theme","Oth_ClassII")
'''
colnames(c2_smm_asso) <- c("RecInfo","DonInfo","Abv_ClassII","Total_eps_ClassII","DRB","AbDRB","otDRB","DQB",
                           "AbDQB","otDQB","DQA","AbDQA","otDQA","DPB","AbDPB","otDPB",
                           "DPA","AbDPA","otDPA","theme","Oth_ClassII")




'
"RecInfo","DonInfo","All_ABV_ClassII","Total_eps","allDRB","AbDRB","otDRB","allDQB",
"AbDQB","otDQB","allDQA","AbDQA","otDQA","allDPB","AbDPB","otDPB",
"allDPA","AbDPA","otDPA","theme","All_Oth_ClassII"
'
c1_asso %>% merge(c2_smm_asso) %>% 
  mutate('Total_eps' = Total_eps_ClassI + Total_eps_ClassII,
         'Total_Abv' = Abv_ClassI + Abv_ClassII,
         'Total_Oth' = Oth_ClassI + Oth_ClassII) -> c1_c2_eplet_forAsso
head(c1_c2_eplet_forAsso)
#c1_c2_eplet_forAsso %>% writexl::write_xlsx("Association/new_pheno_cov/c1c2_eplet_forAsso_withSMMS_1147.xlsx")


c1_asso %>% merge(c2_smm_asso) %>% 
  mutate('Total_eps' = Total_eps_ClassI + Total_eps_ClassII,
         'Total_Abv' = Abv_ClassI + Abv_ClassII,
         'Total_Oth' = Oth_ClassI + Oth_ClassII) %>%
  merge(pheno,by.x = "RecInfo",by.y = "bCODE") -> c1_c2_eplet_forAsso

ncol(c1_c2_eplet_forAsso)
head(c1_c2_eplet_forAsso)
colnames(c1_c2_eplet_forAsso)

#c1_c2_eplet_forAsso %>% writexl::write_xlsx("Association/c1c2_eplet_forAsso_withSMMS.xlsx")
#c1_c2_eplet_forAsso %>% writexl::write_xlsx("Association/new_pheno_cov/c1c2_eplet_forAsso_withSMMS_1147.xlsx")
#  filter(!(theme == "SMMS" & !( %in% c("Total_eps","Total_Abv","Total_Oth","Total_eps_ClassII","Abv_ClassII","Oth_ClassII","AbDQB","AbDRB","DPB","DRB","otDQB","otDRB"))))

c1_c2_eplet_forAsso %>% filter(theme == "EMS") %>% select(-theme)
#c1_c2_eplet_forAsso %>% filter(theme == "SMMS") %>% select(-theme)
#c2_data_scale[,i] <- scale(as.numeric(c2_data_scale[,i]))
result <- NULL
for(i in 3:27){
  if (colnames(c1_c2_eplet_forAsso)[i] != "theme") {
    c1_c2_eplet_forAsso[,i] <- scale(as.numeric(c1_c2_eplet_forAsso[,i]))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme == "EMS"), family="binomial")
    result <- rbind(result, c("EMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme != "EMS"), family="binomial")
    result <- rbind(result, c("SMMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
  }  
}
head(result)

colnames(result) <- c("theme","ID","Estimate","SE","Z","P")

result <- as.data.frame(result)
for (i in 3:6) {
  result[,i] <- as.numeric(result[,i])
 
}

head(result)
str(result)
result$ID
result %>% 
  mutate('log10P' = log10(as.numeric(P))) %>% #head()
  mutate("class" = ifelse(ID %in% c("Total_eps_ClassI","Abv_ClassI","Oth_ClassI"),"CLASS I",
                          ifelse(ID %in% c("Total_eps","Total_Abv","Total_Oth"),"CLASS I + II",
                          ifelse(ID %in% c("Total_eps_ClassII","Abv_ClassII","Oth_ClassII"),"CLASS II","CLASS II (Gene)")))) %>% #head()
  #filter(!(theme == "SMMS" & !(ID %in% c("Total_eps","Total_Abv","Total_Oth","Total_eps_ClassII","Abv_ClassII","Oth_ClassII","AbDPB","AbDRB","allDPB","allDRB","otDPB","otDRB")))) -> result_forplot
filter(!(theme == "SMMS" & !(ID %in% c("Total_eps","Total_Abv","Total_Oth","Total_eps_ClassII","Abv_ClassII","Oth_ClassII","AbDQB","AbDRB","DPB","DQB","otDQB","otDRB")))) -> result_forplot


#writexl::write_xlsx(result_forplot,"Association/c1c2_eplet_asso_withouthHLAms_scale_withSMMS_onlyEpletCount.xlsx")


df <- read_xlsx("Association/c1c2_eplet_asso_withouthHLAms_scale_withSMMS_onlyEpletCount.xlsx")
head(df)
table(df$ID)
df %>% mutate(ID = ifelse(grepl("Total_eps",ID),"ALL",
                   ifelse(grepl("Abv",ID),"Abv",
                   ifelse(grepl("Oth",ID),"Oth",ID)))) -> df

head(df)
  #ggplot(aes(x=ID,y=-log10P,color=factor(class,levels = c("CLASS I + II","CLASS I","CLASS II")),shape=theme,size=Estimate)) + 
ggplot(df,aes(x=ID,y=-log10P,shape=theme,colour=Estimate)) + 
  geom_point(size = 2) +
  scale_color_gradient2(low="blue",mid ="green",high="red",midpoint =0.1) +
  labs(colour='Estimate',
       shape = "")  +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        #legend.title = element_blank(),
        legend.text = element_text(size = 13)) + 
  theme(strip.text.x = element_text(size = 13,face = "bold")) + 
  facet_grid(~factor(class,levels = c("CLASS I + II","CLASS I","CLASS II","CLASS II (Gene)")),scales = "free_x",space = "free_x") + 
  geom_text_repel(data=subset(df,-log10P > 1.25),
                  aes(label=ID),
                  show.legend = F,
                  size =4,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))
  #theme_bw() 



'''
ggplot(c2_result_v2,aes(x=fct_inorder(ID),y=-log10P,color=Gene,shape=type)) +
  geom_point(size = 2) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom") + 
  geom_text_repel(data=subset(c2_result_v2,-log10P > 1.25),
                  aes(label=ID),
                  show.legend = F,
                  size =4,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))
'''
  
