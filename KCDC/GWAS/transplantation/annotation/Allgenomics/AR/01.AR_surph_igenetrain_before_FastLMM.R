#### AR asso for igeneTrain


## phenotype
library(tidyverse)
setwd("~/Desktop/KCDC/transplantation/allogenomic/")

pheno <- readxl::read_xlsx("phenotype/DATASET_20250315_Final.xlsx")
head(pheno)

ref <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/00.sampleInfo/ALL(2019and2020).sampleID.withbCODE.xlsx")
head(ref)
qcin <- readxl::read_xlsx("../00.sampleInfo/KOTRY.KBAv1.1????????????????????????????????????_20250311_v2.xlsx",sheet = 2)
head(qcin)

ref %>% filter(bCODE %in% qcin$bCODE) %>% dim()

final_sample_list <- read.table("2025_AR/KBA_QC/2nd/change.ID.finalQCin.list.txt")

head(final_sample_list)

ref %>% filter(bCODE %in% qcin$bCODE) -> ref

final_sample_list %>% filter(V2 %in% ref$KBA_ID) %>% dim()

head(ref)
ref %>% filter(type == "KR") %>% select(-Prod,-type) -> KR

ref %>% filter(type == "KD") %>% select(-Prod,-type,-OriID) -> KD

head(KR)
head(KD)
colnames(KR) <- c("KBA_ID.KR","bCODE.KR","OriID","ref")
colnames(KD) <- c("KBA_ID.KD","bCODE.KD","ref")

head(KR)
head(KD)
head(final_sample_list)
KR %>% left_join(KD) %>% filter(KBA_ID.KR %in% final_sample_list$V2) %>%
  filter(KBA_ID.KD %in% final_sample_list$V2) -> out
#1441
head(pheno)
out %>% count(OriID %in% pheno$SUBJNO)
head(out)

writexl::write_xlsx(out,"2025_AR/KBA_QC/KBA.QC.list.AR_1441pairs.20250325.xlsx")
pheno %>% filter(SUBJNO %in% out$OriID) %>% 
  mutate(AR_BX = ifelse(is.na(AR_BX_FIRST_DATE),0,1)) %>%
  mutate(GF = ifelse(is.na(GF_DATE),0,1)) %>%
  count(is.na(AR_BX_FIRST_DATE))

pheno %>% filter(SUBJNO %in% out$OriID) %>% count(is.na(AR_BX_FIRST_DATE))
pheno %>% filter(SUBJNO %in% out$OriID) %>% count(is.na(AR_TREAT_FIRST_DATE))
pheno %>% filter(SUBJNO %in% out$OriID) %>% count(is.na(GF_DATE))
pheno %>% filter(SUBJNO %in% out$OriID) %>% count(is.na(GF_MODALITY))
pheno %>% filter(SUBJNO %in% out$OriID) %>% head()

################################################################################
### survival analysis


library(survival)
library(survminer)  # 예쁜 그래프!
library(tidyverse)
setwd("~/Desktop/KCDC/transplantation/allogenomic/")



AR_pheno <- readxl::read_xlsx("phenotype/DATASET_20250315_Final.xlsx")
pca <- read.table("2025_AR/iGeneTrain_Asso/PCA.txt",header = T)
head(AR_pheno)
head(pca)
ref <- readxl::read_xlsx("../00.sampleInfo/ALL(2019and2020).sampleID.withbCODE.xlsx")
head(pca)
head(AR_pheno)
head(ref)

pca %>% left_join(ref %>% select(KBA_ID,bCODE) %>% rename(IID = KBA_ID)) -> pca


AR_pheno %>%  mutate(AR_TREAT = ifelse(is.na(AR_TREAT_FIRST_DATE),0,1)) %>% #select(KT_DATE,AR_TREAT_FIRST_DATE,LAST_FU_DATE)
  mutate(AR_BX = ifelse(is.na(AR_BX_FIRST_DATE),0,1)) %>% #select(KT_DATE,AR_TREAT_FIRST_DATE,LAST_FU_DATE)
  mutate(GF = ifelse(is.na(GF_DATE),0,1)) %>% #select(KT_DATE,AR_TREAT_FIRST_DATE,LAST_FU_DATE)
  mutate(AR = ifelse(AR_TREAT==1 | AR_BX ==1,1,0)) %>% #count(AR,AR_TREAT,AR_BX) 
  mutate(KT_DATE = as.Date(KT_DATE),
         AR_TREAT_FIRST_DATE = as.Date(AR_TREAT_FIRST_DATE),
         AR_BX_FIRST_DATE = as.Date(AR_BX_FIRST_DATE),
         GF_DATE = as.Date(GF_DATE),
         LAST_FU_DATE = as.Date(LAST_FU_DATE)) %>% #select(KT_DATE,AR_TREAT_FIRST_DATE,LAST_FU_DATE)
  mutate(time_event_day_AR_TREAT  = case_when(
    !is.na(AR_TREAT_FIRST_DATE) ~ as.numeric(AR_TREAT_FIRST_DATE - KT_DATE, units = "days"),
    TRUE ~ as.numeric(LAST_FU_DATE - KT_DATE, units = "days"))) %>%
  mutate(time_event_day_AR_BX  = case_when(
    !is.na(AR_BX_FIRST_DATE) ~ as.numeric(AR_BX_FIRST_DATE - KT_DATE, units = "days"),
    TRUE ~ as.numeric(LAST_FU_DATE - KT_DATE, units = "days"))) %>% 
  mutate(time_event_day_GF  = case_when(
    !is.na(GF_DATE) ~ as.numeric(GF_DATE - KT_DATE, units = "days"),
    TRUE ~ as.numeric(LAST_FU_DATE - KT_DATE, units = "days"))) -> pheno

pca %>% left_join(pheno %>% rename(bCODE = ...1)) -> final_pheno


final_pheno %>% 
  summarise(across(everything(), ~ n_distinct(.))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_levels") %>%
  filter(n_levels == 1)


#AR_BX_FIRST_DATE
#AR_TREAT_FIRST_DATE
#GF_DATE
# age, gender, donor age, donor gender, living/deceased donor, year of transplant 
# age, gender, donor age, donor gender, year of transplant, living/deceased donor, and acute rejection 
#table(final_pheno$SUBJNO)


factor_vars <- c("SEX", "D_SEX", "DONOR_TYPE")
final_pheno <- final_pheno %>%
  mutate(across(all_of(factor_vars), as.factor))


colnames(final_pheno)
coxph_model_AR_TREAT <- coxph(
  Surv(time_event_day_AR_TREAT, AR_TREAT) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO,
  data = final_pheno
)
#donor type 다 living


coxph_model_AR_BX <- coxph(
  Surv(time_event_day_AR_BX, AR_BX) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO,
  data = final_pheno
)

coxph_model_GF <- coxph(
  Surv(time_event_day_GF, GF) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO + AR,
  data = final_pheno
)



residuals_AR_TREAT <- residuals(coxph_model_AR_TREAT, type = "martingale")
residuals_AR_BX <- residuals(coxph_model_AR_BX, type = "martingale")
residuals_GF <- residuals(coxph_model_GF, type = "martingale")

residual_df <- data.frame(IID = final_pheno$IID,AR_TREAT = residuals_AR_TREAT,AR_BX = residuals_AR_BX,GF = residuals_GF)
residual_df %>% head()

write.table(residual_df,"2025_AR/iGeneTrain_Asso/coxph.residuals.forAR.txt",col.names = T,row.names = F,quote = F,sep = "\t")

## surplot
# GF를 그룹 변수로
final_pheno$GF_group <- factor(final_pheno$GF, labels = c("No GF", "GF"))

# 생존 객체
surv_object <- Surv(final_pheno$time_event_day_AR_TREAT, final_pheno$AR_TREAT)

# 그룹에 따라 생존곡선 적합
fit <- survfit(surv_object ~ GF_group, data = final_pheno)

# 시각화
ggsurvplot(
  fit,
  data = final_pheno,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.title = "Graft Failure",
  legend.labs = c("No GF", "GF"),
  palette = "Dark2"
)




