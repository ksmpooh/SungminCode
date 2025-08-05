### suvival
library(tidyverse)
library(ggrepel)
library(scales)
library(survival)
library(survminer)
library(cowplot)

setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/ori_pair/")
#~/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/ori_pair/KRKD.SIRP-a.1148pair.txt
pheno <- readxl::read_xlsx("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso_v2.xlsx")

sample_ref <- read.table("~/Desktop/KCDC/transplantation/allogenomic/02.KR.KD.immune.cell.co-signal_targetGene.alleleMatching01.Score_Sum.txt",header = T)
head(sample_ref) 
sample_ref %>% select(1,2) -> ref

sirp <- read.table("KRKD.SIRP-a.1148pair.txt",header = T)
head(sirp)
sirp %>% select(1,4) -> sirp

library(survival)

#### KR KD 
head(ref)
head(sirp)
colnames(sirp) <- c("KBA_ID.KR","KR.hap")
ref  %>% left_join(sirp) -> df
head(df)
colnames(sirp) <- c("KBA_ID.KD","KD.hap")
df  %>% left_join(sirp) -> ori_df
df <-ori_df
head(df)
df %>% mutate(sirp_missmatch = ifelse(KR.hap == KD.hap,0,1)) -> df
table(df$sirp_missmatch)
head(df)
pheno %>% left_join(df %>% select(2,5)) -> tmp
head(tmp)
tmp %>% select(rej_time_day,rej_tot)
tmp %>% count(rej_tot,sirp_missmatch)

head(tmp)

head(df)
pheno %>% left_join(df) %>% select(RecInfo,DonInfo,KR.hap,KD.hap) %>% rename(SIRPa_type.KR = KR.hap,SIRPa_type.KD = KD.hap) %>% 
  write.table("KRKD.SIRPa_type.1148pair.txt",col.names = T,row.names = F,quote = F,sep = "\t")


coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~factor(sirp_missmatch)+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% as.data.frame())
coxfit
ggadjustedcurves(coxfit, data=tmp %>% as.data.frame(), variable="sirp_missmatch",xlab="Month",legend.title = "SIRP-a")  -> a
a
result <- rbind(result,c("EMS","Total_eps","basic",as.vector(summary(coxfit)$coefficients[1,])))
result

coxfit = survfit(Surv(rej_time_day/30,rej_tot)~sirp_missmatch+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% as.data.frame())

ggsurvplot(coxfit,ggsu)
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), 
           conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "SIRP-a") 
#  xlab("Month")


ori_df #%>% mutate(sirp_missmatch = ifelse(KR.hap == KD.hap,0,1)) %>%
  
ori_df %>%
  # KR.hap -> kr1, kr2로 분리
  separate(KR.hap, into = c("kr1", "kr2"), sep = "/") %>%
  # KD.hap -> kd1, kd2로 분리
  separate(KD.hap, into = c("kd1", "kd2"), sep = "/") %>%
  rowwise() %>%
  mutate(
    sirp_missmatch = 2 - sum(c(kr1, kr2) %in% c(kd1, kd2))
  ) %>%
  ungroup() -> df

head(df)
pheno %>% left_join(df %>% select(2,7)) -> tmp
head(tmp)

coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~sirp_missmatch+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% as.data.frame())
coxfit
ggadjustedcurves(coxfit, data=tmp %>% as.data.frame(), variable="sirp_missmatch",xlab="Month",legend.title = "SIRP-a missmatch")  -> a
a

coxfit

head(pheno)
pheno %>% select(rej_tot,rej_time_day) %>% mutate(yer = rej_time_day/365) %>% mutate(check = ifelse(yer < 1,1,0)) %>%
  count(rej_tot,check)

pheno %>% select(rej_tot,rej_time_day) %>% mutate(month = round(rej_time_day/30,digits = 0)) %>% #head()
  ggplot(aes(x=month,fill=factor(rej_tot))) +
  #geom_bar(stat = 'identity',position = "dodge")
  geom_bar(position = "dodge")

a <- coxph(Surv(rej_time_day,rej_tot)~1,data=pheno)
summary(a)
coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~factor(sirp_missmatch)+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% as.data.frame())
coxfit






######### new_pair 1441 20250423




setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/new_pair/")
#~/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/ori_pair/KRKD.SIRP-a.1148pair.txt
pheno <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/allogenomic/phenotype/DATASET_20250315_Final.xlsx")
sirp <- read.table("KRKD.KOTRY.AR.SIRPA.txt",header = T)
head(sirp)
ref <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/00.sampleInfo/ALL(2019and2020).sampleID.withbCODE.xlsx")

pheno %>%  mutate(AR_TREAT = ifelse(is.na(AR_TREAT_FIRST_DATE),0,1)) %>% #select(KT_DATE,AR_TREAT_FIRST_DATE,LAST_FU_DATE)
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
    TRUE ~ as.numeric(LAST_FU_DATE - KT_DATE, units = "days"))) -> AR_pheno

head(AR_pheno)
head(sirp)
sirp %>% left_join()pheno %>% rename(bCODE = ...1)
head(ref)

pca <- read.table("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/iGeneTrain_Asso/PCA.txt",header = T)

sirp %>% left_join(ref %>% filter(type == "KR") %>% rename(KBA_ID.KR = KBA_ID) %>% select(KBA_ID.KR,bCODE)) %>% 
  left_join(AR_pheno %>% rename(bCODE = ...1)) %>%
  left_join(pca %>% rename(KBA_ID.KR = FID)) -> final_pheno
head(sirp)


#setwd("~/Desktop/KCDC/transplantation/allogenomic/")

head(final_pheno)


AR_pheno %>% 
  summarise(across(everything(), ~ n_distinct(.))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_levels") %>%
  filter(n_levels == 1)

#############
head(final_pheno)
final_pheno %>% #count(SIRPa_type.KR,SIRPa_type.KD) %>%
  filter(SIRPa_type.KR %in% c("v1/v1","v1/v2","v2/v2")) %>%
  filter(SIRPa_type.KD %in% c("v1/v1","v1/v2","v2/v2")) %>% 
  mutate(SIRPa_path_DtoR = paste0(SIRPa_type.KD," -> ",SIRPa_type.KR)) -> final_pheno_onlyv1v2
  

final_pheno %>% count(SIRPa_type.KR) %>% 
final_pheno %>% count(SIRPa_type.KD)

dim(final_pheno_onlyv1v2)
final_pheno_onlyv1v2 %>% count(SIRPa_type.KR)
final_pheno_onlyv1v2 %>% count(SIRPa_type.KD)



coxph_model_AR_TREAT
summary(coxph_model_AR_TREAT)
coxph_model_AR_TREAT <- coxph(
  Surv(time_event_day_AR_TREAT, AR_TREAT) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO,
  data = final_pheno_onlyv1v2
)
#donor type 다 living


coxph_model_AR_BX <- coxph(
  Surv(time_event_day_AR_BX, AR_BX) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO,
  data = final_pheno_onlyv1v2
)

coxph_model_GF <- coxph(
  Surv(time_event_day_GF, GF) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO + AR,
  data = final_pheno_onlyv1v2
)

library(survival)
library(survminer)
surv_object <- Surv(final_pheno_onlyv1v2$time_event_day_GF, final_pheno_onlyv1v2$GF)
#Head(final_pheno)
# 그룹 기준 생존곡선 추정
fit <- survfit(surv_object ~ SIRPa_path_DtoR, data = final_pheno_onlyv1v2)

# 생존곡선 시각화
ggsurvplot(
  fit,
  data = final_pheno_onlyv1v2,
  pval = TRUE,                  # p-value 표시
  conf.int = TRUE,    신뢰구간
  risk.table = TRUE,            # 생존자 수 테이블
  legend.title = "SIRPα D→R",   # 범례 제목
  legend.labs = levels(as.factor(final_pheno_onlyv1v2$SIRPa_path_DtoR)),
  palette = "Dark2"             # 색상 팔레트
)


library(dplyr)
library(tidyr)
library(rlang)

make_sirpa_table <- function(data, phenotype_col_name) {
  # 동적 변수 심볼화
  pheno <- sym(phenotype_col_name)
  
  # 1?????? 비율 계산
  summary <- data %>%
    group_by(SIRPa_type.KD, SIRPa_type.KR) %>%
    summarise(
      count_1 = sum(!!pheno == 1, na.rm = TRUE),
      total = sum(!is.na(!!pheno)),
      .groups = "drop"
    ) %>%
    mutate(ratio_str = paste0(count_1, "/", total, " (", round(count_1 / total * 100, 1), "%)"))
  
  # 2?????? 3x3 테이블로 변환
  ratio_table <- summary %>%
    select(SIRPa_type.KD, SIRPa_type.KR, ratio_str) %>%
    pivot_wider(names_from = SIRPa_type.KR, values_from = ratio_str) %>%
    rename(`Donor (KD)` = SIRPa_type.KD)
  
  # 3?????? 행 Total
  row_total <- summary %>%
    group_by(SIRPa_type.KD) %>%
    summarise(
      count_1 = sum(count_1),
      total = sum(total),
      .groups = "drop"
    ) %>%
    mutate(Total = paste0(count_1, "/", total, " (", round(count_1 / total * 100, 1), "%)")) %>%
    select(SIRPa_type.KD, Total)
  
  ratio_table <- left_join(ratio_table, row_total, by = c("Donor (KD)" = "SIRPa_type.KD"))
  
  # 4?????? 열 Total
  col_total <- summary %>%
    group_by(SIRPa_type.KR) %>%
    summarise(count_1 = sum(count_1), total = sum(total), .groups = "drop") %>%
    mutate(ratio_str = paste0(count_1, "/", total, " (", round(count_1 / total * 100, 1), "%)")) %>%
    pivot_wider(names_from = SIRPa_type.KR, values_from = ratio_str)
  
  # 5?????? 전체 Total
  overall <- summary %>%
    summarise(count_1 = sum(count_1), total = sum(total)) %>%
    mutate(overall_ratio = paste0(count_1, "/", total, " (", round(count_1 / total * 100, 1), "%)")) %>%
    pull(overall_ratio)
  
  # 6?????? Total 행 생성 (모든 열을 character로)
  total_row <- c("Total", col_total, overall) %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(across(everything(), as.character))
  
  # 7?????? 결합 및 반환
  final_table <- bind_rows(ratio_table, total_row)
  
  return(final_table)
}
final_pheno_onlyv1v2
make_sirpa_table(final_pheno_onlyv1v2, "AR_TREAT") %>% writexl::write_xlsx("SIRP.direction.AR_TRAET.count.xlsx")
make_sirpa_table(final_pheno_onlyv1v2, "AR_BX") %>% writexl::write_xlsx("SIRP.direction.AR_BX.count.xlsx")
make_sirpa_table(final_pheno_onlyv1v2, "GF") %>% writexl::write_xlsx("SIRP.direction.GF.count.xlsx")


head(final_pheno_onlyv1v2)
final_pheno %>% count(AR_TREAT)
final_pheno %>% count(AR_BX)
final_pheno %>% count(GF)

#########
final_pheno
head(final_pheno)

coxph_model_GF <- coxph(
  Surv(time_event_day_GF, GF) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO + AR,
  data = final_pheno
)

final_pheno
# Cox 모델에 따라 생존곡선 적합
surv_fit_GF <- survfit(
  Surv(time_event_day_GF, GF) ~ sirp_missmatch,
  data = final_pheno
)

ggsurvplot(
  surv_fit_GF,
  data = final_pheno,
  risk.table = TRUE,          # 하단에 위험표 표시
  pval = TRUE,                # p-value 표시
  conf.int = TRUE,            # 신뢰구간 표시
  palette = "Dark2",          # 색상 팔레트
  legend.title = "SIRPα mismatch",
  legend.labs = c("Match", "Mismatch"),  # 범례 라벨
  xlab = "Days since transplant",
  ylab = "Graft survival probability",
  title = "Graft Survival by SIRPα mismatch",
  font.title = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.legend = c(12)
)

# 생존곡선 시각화
ggsurvplot(
  fit,
  data = final_pheno_onlyv1v2,
  pval = TRUE,                  # p-value 표시
  conf.int = TRUE,              # 신뢰구간
  risk.table = TRUE,            # 생존자 수 테이블
  legend.title = "SIRPα D→R",   # 범례 제목
  legend.labs = levels(as.factor(final_pheno_onlyv1v2$SIRPa_path_DtoR)),
  palette = "Dark2"             # 색상 팔레트
)


cox_model <- coxph(Surv(time_event_day_GF, GF) ~ AGE + SEX + sirp_missmatch, data = final_pheno)

# 그룹별 생존곡선 예측용 데이터
new_data <- data.frame(
  AGE = mean(final_pheno$AGE),
  SEX = "Male",
  sirp_missmatch = factor(c("Match", "Mismatch"), levels = c("Match", "Mismatch"))
)

# 예측
surv_fit <- survfit(cox_model, newdata = new_data)

# 시각화
ggsurvplot(surv_fit, data = new_data, legend.labs = c("Match", "Mismatch"))


ggsurvplot(
  surv_fit,
  data = final_pheno,
  risk.table = TRUE,          # 하단에 위험표 표시
  pval = TRUE,                # p-value 표시
  conf.int = TRUE,            # 신뢰구간 표시
  palette = "Dark2",          # 색상 팔레트
  legend.title = "SIRPα mismatch",
  legend.labs = c("Match", "Mismatch"),  # 범례 라벨
  xlab = "Days since transplant",
  ylab = "Graft survival probability",
  title = "Graft Survival by SIRPα mismatch",
  font.title = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.legend = c(12)
)

# 1. factor 정리
final_pheno
# 2. 모델 생성
cox_model <- coxph(Surv(time_event_day_GF, GF) ~ AGE + SEX + , data = final_pheno)
final_pheno$sirp_missmatch <- factor(final_pheno$sirp_missmatch, levels = c("Match", "Mismatch"))

# 3. 예측 대상
new_data <- data.frame(
  AGE = mean(final_pheno$AGE, na.rm = TRUE),
  #SEX = "Male",
  sirp_missmatch = factor(c("Match", "Mismatch"), levels = levels(final_pheno$sirp_missmatch))
)

# 4. 생존곡선 예측
surv_fit <- survfit(cox_model, newdata = new_data)

# 5. 시각화
ggsurvplot(
  surv_fit,
  data = new_data,  # <- 여기서 `final_pheno`가 아닌 `new_data`
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  palette = "Dark2",
  legend.title = "SIRPα mismatch",
  legend.labs = c("Match", "Mismatch"),
  xlab = "Days since transplant",
  ylab = "Graft survival probability",
  title = "Graft Survival by SIRPα mismatch",
  font.title = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.legend = c(12)
)
pheno
head(final_pheno)
final_pheno %>% rename(bCODE.KR = bCODE,bCODE.KD=...2) %>%
  select(bCODE.KR,bCODE.KD,SIRPa_type.KR,SIRPa_type.KD,sirp_missmatch,PC1:PC10) %>% writexl::write_xlsx("KOTRY.AR.1441pair.SIRP_a_haplotype.withPC.xlsx")


head(final_pheno)
final_pheno


#load("/Users/ksmpooh/Downloads/DATASET20250425.rds")
data_object <- readRDS("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/phenotype/DATASET20250425.rds")
table(data_object$DONOR_TYPE)

data_object
