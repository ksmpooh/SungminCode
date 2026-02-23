#### PIRCHE vs HLA MM
library(tidyverse)
library(stringr)
library(readxl)

setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/")

pheno <- read_xlsx("~/Desktop/KCDC/transplantation/allogenomic/phenotype/DATASET_20250315_Final.xlsx") %>%
  mutate(AR_TREAT = ifelse(is.na(AR_TREAT_FIRST_DATE),0,1)) %>% #select(KT_DATE,AR_TREAT_FIRST_DATE,LAST_FU_DATE)
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
    TRUE ~ as.numeric(LAST_FU_DATE - KT_DATE, units = "days")))

pirche <- read_xlsx("PIRCHE/pirche_result_KOTRY_250421.xlsx",skip = 1)
head(pirche)
colnames(pirche)

pirche %>% filter(!is.na(Patient_Donor_ID)) %>% filter(!is.na(PIRCHE_I)) %>% 
  select(Patient_Donor_ID,PIRCHE_II,A_Originated_Epitopes_count:DRB1_Originated_Epitopes_count,DQB1_Originated_Epitopes_count,DPB1_Originated_Epitopes_count) %>%
  pivot_longer(PIRCHE_II:DPB1_Originated_Epitopes_count) %>% 
  ggplot(aes(x=name,y=value)) + 
  geom_violin()


pirche %>% filter(!is.na(Patient_Donor_ID)) %>% filter(!is.na(PIRCHE_I)) %>% 
  select(Patient_Donor_ID,PIRCHE_II,A_Originated_Epitopes_count:DRB1_Originated_Epitopes_count,DQB1_Originated_Epitopes_count,DPB1_Originated_Epitopes_count) %>%
  mutate(
    PIRCHE_classI = rowSums(across(c(A_Originated_Epitopes_count,
                                     B_Originated_Epitopes_count,
                                     C_Originated_Epitopes_count)),
                            na.rm = TRUE),
    PIRCHE_classII = rowSums(across(c(DRB1_Originated_Epitopes_count,
                                      DQB1_Originated_Epitopes_count,
                                      DPB1_Originated_Epitopes_count)),
                             na.rm = TRUE)) -> pirche_pro
  #count(PIRCHE_II == PIRCHE_classI + PIRCHE_classII ) 
  #select(Patient_Donor_ID,PIRCHE_II,PIRCHE_classI,PIRCHE_classII) -> 
  #filter(PIRCHE_II != PIRCHE_classI + PIRCHE_classII )

head(pheno)
pirche_pro  %>% mutate(SUBJNO = str_replace_all(Patient_Donor_ID,"dnr0","KR")) %>% 
  left_join(pheno)
pirche_pro

c1 <- read_excel("hla_matchmaker/c1_forWAS.xlsx")
c2 <- read_excel("hla_matchmaker/c2_forWAS.xlsx")
head(c1)
head(c2)
colnames(c1)
colnames(c2)


c1 %>% select(RecInfo,`All Eplets`,AbVEp,OthEp) %>% rename(HLAmm_classI=`All Eplets`) %>%
  left_join(c2 %>% mutate(All_Other_ClassII = Total_eps - All_ABV_ClassII) %>% 
              select(RecInfo,All_ABV_ClassII,Total_eps,All_Other_ClassII) %>%
              rename(HLAmm_classII = Total_eps)) %>% 
  mutate(HLAmm_all = HLAmm_classI + HLAmm_classII) -> hlamm

head(pirche_pro)
head(pheno)
pirche_pro %>% tail()
pheno %>% select(SUBJNO) %>% tail()
pirche_pro %>%mutate(SUBJNO = str_replace_all(Patient_Donor_ID,"dnr","KR")) %>% select(SUBJNO)

hlamm %>% left_join(pheno %>% rename(RecInfo = bCODE_R) %>% select(RecInfo,SUBJNO)) %>% #head()
  left_join(pirche_pro %>%mutate(SUBJNO = str_replace_all(Patient_Donor_ID,"dnr","KR")) %>% 
              left_join(pheno) %>% select(SUBJNO,PIRCHE_II,PIRCHE_classI,PIRCHE_classII)) %>%
  select(RecInfo,HLAmm_classI,HLAmm_classII,HLAmm_all,PIRCHE_II,PIRCHE_classI,PIRCHE_classII) %>% #head()
  pivot_longer(
    cols = HLAmm_classI:PIRCHE_classII,
    names_to = "name",
    values_to = "score"
  ) %>% #head
  mutate(method = if_else(str_starts(name, "HLAmm"), "HLAmm", "PIRCHE")) %>%
  mutate(type = case_when(
    str_detect(name,"classII") ~ "classII",
    str_detect(name,"classI") ~ "classI",
    TRUE ~"Total")) %>%
  select(RecInfo, method,type, score) -> hlamm_pirche
  # 같은 type 내에서 HLAmm, PIRCHE를 나란히 보게 wide화

head(hlamm_pirche)
head(pheno)
hlamm_pirche %>% 
  left_join(pheno %>% rename(RecInfo = bCODE_R) %>% select(RecInfo,AR_TREAT,AR_BX,AR)) %>% #head
  select(RecInfo,method:score,AR_TREAT,AR_TREAT,AR_BX,AR) %>% #head()
  pivot_longer(
    cols = c(AR_TREAT, AR_BX, AR),
    names_to = "Condition", values_to = "Status"
  ) %>% 
  mutate(
    # 보기 좋은 순서/라벨
    Condition = factor(Condition, levels = c("AR_TREAT","AR_BX","AR")),
    Status    = factor(Status, levels = c(0,1), labels = c("Control","Case"))
  ) %>% #head()
  ggplot(aes(x=type,y=score,fill=Status,color=method)) + 
  geom_violin() + 
    facet_grid(~Condition, scales = "free") +
    scale_fill_manual(values = c("#bdbdbd", "#3182bd"), name = "Status") +
    labs(
      x = "Score", y = "Sample count",
      title = "Distributions of Eplet metrics by AR conditions (Class II)"
    ) +
    theme_bw()

#1441 * 2 * 3 * 3
temp
hlamm_pirche %>%
  pivot_wider(
    names_from = method,
    values_from = score
  ) %>% head()
  ggplot(aes(x=HLAmm,y=PIRCHE)) + 
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ type, scales = "free") +
  labs(x = "HLA Matchmaker (eplet count / all, classI, classII)",
       y = "PIRCHE score (II, classI, classII)",
       title = "HLAmm vs PIRCHE by Pairing Type") + 
  theme_bw()

çççc



   
 library(dplyr)
library(ggplot2)

# 1️⃣ 각 type별 상관계수 계산
cor_df <- hlamm_pirche %>%
  pivot_wider(names_from = method, values_from = score) %>%
  group_by(type) %>%
  summarise(
    cor_r = cor(HLAmm, PIRCHE, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

# 2️⃣ 원 데이터
plot_df <- hlamm_pirche %>%
  pivot_wider(names_from = method, values_from = score)

# 3️⃣ 시각화
ggplot(plot_df, aes(x = HLAmm, y = PIRCHE)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  #coord_fixed(ratio = 1) + 
  facet_wrap(~ type, scales = "free") + 
  # annotate로 상관계수 추가
  geom_text(
    data = cor_df,
    aes(x = Inf, y = Inf, 
        label = paste0("Spearman r = ", round(cor_r, 3))),
    hjust = 1.1, vjust = 1.3,
    size = 3.5,
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = "HLA Matchmaker (eplet count / all, classI, classII)",
    y = "PIRCHE score (II, classI, classII)",
    title = "HLAmm vs PIRCHE by Pairing Type"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

### norm

hlamm_pirche %>% group_by(method,type) %>%
  mutate(norm_score = (score - min(score, na.rm = TRUE)) / 
           (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))) %>%
  ungroup() %>% select(-score) %>%
  pivot_wider(
    names_from = method,
    values_from = norm_score,
    ) %>%
  ggplot(aes(x=HLAmm,y=PIRCHE)) + 
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ type, scales = "free") +
  labs(x = "HLA Matchmaker (eplet count / all, classI, classII)",
       y = "PIRCHE score",
       title = "HLAmm vs PIRCHE by Pairing Type") + 
  theme_bw()


hlamm_norm <- hlamm_pirche %>%
  group_by(method, type) %>%
  mutate(
    .min = min(score, na.rm = TRUE),
    .max = max(score, na.rm = TRUE),
    norm_score = if_else(.max > .min, (score - .min)/(.max - .min), 0)
  ) %>%
  ungroup() %>%
  select(-.min, -.max)

# 2) wide 변환: 정규화된 값으로 HLAmm, PIRCHE 컬럼 만들기
plot_df <- hlamm_norm %>%
  select(-score) %>%
  pivot_wider(
    names_from = method,
    values_from = norm_score
  )

# 3) 각 type별 Spearman 상관계수 계산 (정규화된 값으로)
cor_df <- plot_df %>%
  group_by(type) %>%
  summarise(
    spearman_r = cor(HLAmm, PIRCHE, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

# 4) 시각화: 산점도 + x=y 기준선 + 추세선 + 상관계수
ggplot(plot_df, aes(x = HLAmm, y = PIRCHE)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ type, scales = "free") +
  geom_text(
    data = cor_df,
    aes(x = Inf, y = Inf, label = paste0("Spearman r = ", round(spearman_r, 3))),
    hjust = 1.1, vjust = 1.3, size = 3.5, inherit.aes = FALSE
  ) +
  labs(
    x = "HLA Matchmaker (normalized, all/classI/classII)",
    y = "PIRCHE score (normalized)",
    title = "HLAmm vs PIRCHE by Pairing Type (normalized)"
  ) +
  theme_bw()


####
head(hlamm)
hlamm %>% left_join(pheno %>% rename(RecInfo = bCODE_R) %>% select(RecInfo,SUBJNO)) %>% #head()
  left_join(pirche_pro %>%mutate(SUBJNO = str_replace_all(Patient_Donor_ID,"dnr","KR")) %>% 
              left_join(pheno) %>% select(SUBJNO,PIRCHE_II,PIRCHE_classI,PIRCHE_classII)) ->a
colnames(a) <- c("KR","HLAmm_ClassI_Total","HLAmm_ClassI_AbVEp","HLAmm_ClassI_OthEp","HLAmm_ClassII_AbVEp","HLAmm_ClassII_Total","HLAmm_ClassII_OthEp","HLAmm_Total","sub","PIRCHE_II_Total","PIRCHE_II_classI","PIRCHE_II_classII")

a %>% select(KR,HLAmm_ClassI_Total,HLAmm_ClassI_AbVEp,HLAmm_ClassI_OthEp,HLAmm_ClassII_Total,HLAmm_ClassII_AbVEp,HLAmm_ClassII_OthEp,HLAmm_Total,PIRCHE_II_Total,PIRCHE_II_classI,PIRCHE_II_classII) %>%
  write.xlsx("HLAmm_pirche_KOTRY_MainCategory.2512002_pro.xlsx")
a


pirche_pro %>%mutate(SUBJNO = str_replace_all(Patient_Donor_ID,"dnr","KR")) %>% 
  left_join(pheno) %>% #head()
  select(bCODE_R,2:8) %>% rename(KR = bCODE_R) %>%
  write.xlsx("PIRCHE/pirche_result_KOTRY_2512002_pro.xlsx")


####### association

library(tidyverse)
setwd("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/")

head(pheno)

pirche <- readxl::read_xlsx("PIRCHE/pirche_result_KOTRY_2512002_pro.xlsx")
pirche
main_result <-  readxl::read_xlsx("HLAmm_pirche_KOTRY_MainCategory.2512002_pro.xlsx")

c1 <-read_xlsx("hla_matchmaker/HLAmatchmaker_classI_prop_20251103.xlsx")
c2 <-read_xlsx("hla_matchmaker/HLAmatchmaker_classII_prop_20251103.xlsx")
c2

head(main_result)
head(pheno)
table(pheno$AR_TREAT)
main_result %>% left_join(pheno %>% rename(KR = bCODE_R)) -> main_pheno
colnames(main_result)
main_pheno
#main_pheno
result <- NULL
for(i in 2:ncol(main_result)){
  print(i)
  main_pheno[,i] <- scale(main_pheno[,i])
  temp <- glm(paste("AR_TREAT ~ ",colnames(main_result)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO",sep=""), data=main_pheno, family="binomial")
  result <- rbind(result, c("AR_TREAT",colnames(main_result)[i],as.vector(summary(temp)$coefficients[2,])))
  
  temp <- glm(paste("AR_BX ~ ",colnames(main_result)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO",sep=""), data=main_pheno, family="binomial")
  result <- rbind(result, c("AR_BX",colnames(main_result)[i],as.vector(summary(temp)$coefficients[2,])))
  
  temp <- glm(paste("GF ~ ",colnames(main_result)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO + AR",sep=""), data=main_pheno, family="binomial")
  result <- rbind(result, c("GF",colnames(main_result)[i],as.vector(summary(temp)$coefficients[2,])))
  
  #temp <- glm(paste("rej_tot ~ ",colnames(main_result)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=main_pheno, family="binomial")
  #result <- rbind(result, c("SMMS",colnames(main_result)[i],as.vector(summary(temp)$coefficients[2,])))
}
head(result)
as.vector(summary(temp)$coefficients[2,])
summary(temp)

colnames(result) <- c("phenotype","ID","Estimate","SE","Z","P")

result <- as.data.frame(result)
for (i in 3:6) {
  result[,i] <- as.numeric(result[,i])
  
}
main_result <- result


pirche
#pirche %>% select(1,3:8)
pirche %>% left_join(pheno %>% rename(KR = bCODE_R)) -> pirche_pheno
pirche_pheno

result <- NULL

for(i in 3:ncol(pirche)){
  print(i)
  pirche_pheno[,i] <- scale(pirche_pheno[,i])
  temp <- glm(paste("AR_TREAT ~ ",colnames(pirche_pheno)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO",sep=""), data=pirche_pheno, family="binomial")
  result <- rbind(result, c("AR_TREAT",colnames(pirche_pheno)[i],as.vector(summary(temp)$coefficients[2,])))
  
  temp <- glm(paste("AR_BX ~ ",colnames(pirche_pheno)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO",sep=""), data=pirche_pheno, family="binomial")
  result <- rbind(result, c("AR_BX",colnames(pirche_pheno)[i],as.vector(summary(temp)$coefficients[2,])))
  
  temp <- glm(paste("GF ~ ",colnames(pirche_pheno)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO + AR",sep=""), data=pirche_pheno, family="binomial")
  result <- rbind(result, c("GF",colnames(pirche_pheno)[i],as.vector(summary(temp)$coefficients[2,])))
}

as.vector(summary(temp)$coefficients[2,])
summary(temp)

colnames(result) <- c("phenotype","ID","Estimate","SE","Z","P")
result
result <- as.data.frame(result)
for (i in 3:6) {
  result[,i] <- as.numeric(result[,i])
  
}
pirche_result <- result

head(result)

pirche

head(main_pheno)
###
head(c1)
colnames(c1)
c1 %>% select(RecInfo,AbVEp,OthEp)
colnames(c2)

c1 %>% select(RecInfo,AbVEp,OthEp) %>% left_join(c2 %>% select(RecInfo,Nr_allDRB:Nr_otDPA))

c1 %>% select(RecInfo,AbVEp,OthEp,`1C`:`275K`)%>% rename_with(.cols = 2:ncol(.), .fn = ~ paste0("ClassI_", .x)) %>%
  left_join(c2 %>% select(RecInfo,Nr_allDRB:Nr_otDPA,AbDR_4R:OtPA_190A)) %>% rename(KR = RecInfo) %>% 
  #select(1, where(~ !all(.x == 0)) & 2:391) 
  select(1,2:391 & where(~ !(all(.x == 0) | all(.x == 1)))) -> hlamm

colnames(hlamm)
hlamm %>% select(241:242)
hlamm
hlamm %>% left_join(pheno %>% rename(KR = bCODE_R)) -> eplet_pheno

colnames(eplet_pheno)[i]
eplet_pheno %>% count(`1C`)
eplet_pheno[,i]

result <- NULL
eplet_pheno[, 2:ncol(hlamm)] <- lapply(eplet_pheno[, 2:ncol(hlamm)], function(x) {
  as.numeric(x)})

eplet_pheno[,77]
table(eplet_pheno$classI_44RM)

for(i in 2:ncol(hlamm)){
  #target <- colnames(eplet_pheno)[i]
  print(i)
  eplet_pheno[,i] <- scale(eplet_pheno[,i])
  temp <- glm(paste("AR_TREAT ~ ",colnames(eplet_pheno)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO",sep=""), data=eplet_pheno, family="binomial")
  result <- rbind(result, c("AR_TREAT",colnames(eplet_pheno)[i],as.vector(summary(temp)$coefficients[2,])))
  
  temp <- glm(paste("AR_BX ~ ",colnames(eplet_pheno)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO",sep=""), data=eplet_pheno, family="binomial")
  result <- rbind(result, c("AR_BX",colnames(eplet_pheno)[i],as.vector(summary(temp)$coefficients[2,])))
  
  temp <- glm(paste("GF ~ ",colnames(eplet_pheno)[i], " + AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT + HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH + TRANSPLANT_NO + AR",sep=""), data=eplet_pheno, family="binomial")
  result <- rbind(result, c("GF",colnames(eplet_pheno)[i],as.vector(summary(temp)$coefficients[2,])))
}

as.vector(summary(temp)$coefficients[2,])
summary(temp)

colnames(result) <- c("phenotype","ID","Estimate","SE","Z","P")
result
result <- as.data.frame(result)
for (i in 3:6) {
  result[,i] <- as.numeric(result[,i])
  
}
eplet_result <- result
eplet_result

eplet_result$theme <- "eplet"
pirche_result$theme <- "PIRCHE"
main_result$theme <- "Total"

eplet_result

main_result %>% rbind(pirche_result) %>% rbind(eplet_result) -> merge_df


merge_df
merge_df %>% mutate(Method =
                      case_when(
                        str_detect(ID,"HLAmm") ~ "HLAmatchmaker",
                        str_detect(ID,"PIRCHE") ~ "PIRCHE II",
                        theme == "eplet" ~ "HLAmatchmaker",
                        TRUE ~ "PIRCHE II"          # default
                      )) %>% 
  mutate(type = case_when(
    #str_detect(ID,"HLAmm_Total") ~ str_replace_all(ID,"HLAmm_",""),
    #str_detect(ID,"HLAmm") ~ str_split_fixed(ID,"_",3)[,3],
    str_detect(ID,"HLAmm") ~ str_replace_all(ID,"HLAmm_",""),
    str_detect(ID,"PIRCHE") ~ str_replace_all(ID,"PIRCHE_II_",""),
    str_detect(ID,"count") ~ str_replace_all(ID,"_Originated_Epitopes_count",""),
    str_detect(ID,"classII") ~ str_replace_all(ID,"classII_",""),
    TRUE ~ str_replace_all(ID,"classI_","")          # default
  )) %>% #count(Method,type)
  mutate(class = case_when(
    str_detect(ID,"ClassII") ~ "Class II",
    str_detect(ID,"ClassI") ~ "Class I",
    str_detect(ID,"classII") ~ "Class II",
    str_detect(ID,"classI") ~ "Class I",
    str_detect(ID,"Total") ~ "Total",
    type %in% c("A","B","C") ~ "Class I",
    type %in% c("DRB1","DPA1","DPB1","DQA1","DQB1") ~ "Class II",
    TRUE ~ "Class II")) -> df
head(df)



df %>% filter(theme == 'PIRCHE')
df %>% filter(theme == 'Total')

#df %>% filter(theme == "Total") %>% filter(P < 0.05) %>% select(phenotype,ID,Estimate,SE,Z,P) ^>^ 


df %>% 
  filter(theme == "Total") %>% 
  ggplot(aes(x = type, y = -log10(P),color = phenotype,size = abs(Estimate))) + 
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05),color = "red",linetype = "dashed",linewidth = 0.8) +
  facet_wrap(~ Method, ncol = 2,scales = "free_x", space = "free_x") +
    theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),legend.position = "right") +
  labs(x = "ID",y = "-log10(P)",color = "Phenotype",size = "Effect (|Estimate|)")

head(df)
table(df$theme)

df %>% mutate(Gene = ifelse(theme == "PIRCHE","Gene",ifelse(str_detect(type,"Nr"),"Gene","no"))) %>% filter(Gene == "Gene")
head(df)

df %>% mutate(theme = ifelse(theme == "PIRCHE","Gene",ifelse(str_detect(type,"Nr"),"Gene",theme))) %>% 
  mutate(theme = str_replace_all(theme,"Nr_","")) %>% write_xlsx("AR.associationResult.HLAmm_PIRCHE.xlsx")


df %>% mutate(Gene = ifelse(theme == "PIRCHE","Gene",ifelse(str_detect(type,"Nr"),"Gene","no"))) %>% filter(Gene == "Gene") %>%
  mutate(type = str_replace_all(type,"Nr_","")) %>% filter(P < 0.05) %>% select(phenotype,ID,Estimate,SE,Z,P)

df %>% mutate(Gene = ifelse(theme == "PIRCHE","Gene",ifelse(str_detect(type,"Nr"),"Gene","no"))) %>% 
  mutate(type = str_replace_all(type,"Nr_","")) %>%
  filter(Gene == "Gene") %>%
  ggplot(aes(x = type, y = -log10(P),color = phenotype,size = abs(Estimate))) + 
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05),color = "red",linetype = "dashed",linewidth = 0.8) +
  facet_wrap(~ Method, ncol = 2,scales = "free_x", space = "free_x") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),legend.position = "right") +
  labs(x = "ID",y = "-log10(P)",color = "Phenotype",size = "Effect (|Estimate|)")

df
df %>% filter(theme == "eplet") %>% view
df %>% filter(theme == "eplet") %>% count(class)

df %>% mutate(Gene = ifelse(theme == "PIRCHE","Gene",ifelse(str_detect(type,"Nr"),"Gene","no"))) %>% 
  mutate(type = str_replace_all(type,"Nr_","")) %>%
  filter(Gene != "Gene") %>% filter(theme == "eplet") %>% #count(class)
  mutate(type = str_replace_all(type,"ClassI_","")) %>% filter(P < 0.005) %>% select(phenotype,ID,Estimate,SE,Z,P)# %>% view
  

df %>% mutate(Gene = ifelse(theme == "PIRCHE","Gene",ifelse(str_detect(type,"Nr"),"Gene","no"))) %>% 
  mutate(type = str_replace_all(type,"Nr_","")) %>%
  filter(Gene != "Gene") %>% filter(theme == "eplet") %>% #count(class)
  mutate(type = str_replace_all(type,"ClassI_","")) %>%
  ggplot(aes(x = type, y = -log10(P),color = phenotype,size = abs(Estimate))) + 
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05),color = "red",linetype = "dashed",linewidth = 0.8) +
  geom_text(data = ~ subset(.x, P <= 0.01),aes(label = type),vjust = -0.8,size = 3, color = "black"
  ) +
  facet_wrap(~ class, ncol = 2,scales = "free_x", space = "free_x") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),legend.position = "right",axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  labs(x = "ID",y = "-log10(P)",color = "Phenotype",size = "Effect (|Estimate|)")





coxph_model_AR_TREAT <- coxph(
  Surv(time_event_day_AR_TREAT, AR_TREAT) ~ 
    AGE + SEX + D_AGE + D_SEX + YEAR_OF_TRANSPLANT +
    HLA_A_MISMATCH + HLA_B_MISMATCH + HLA_DR_MISMATCH +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    TRANSPLANT_NO,
  data = final_pheno
)
#donor type �� living
head(pheno)
main_pheno %>% count(AR_TREAT,AR_BX,GF)
main_pheno %>% count(AR_TREAT)
main_pheno %>% count(AR_BX)
main_pheno %>% count(GF)

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
