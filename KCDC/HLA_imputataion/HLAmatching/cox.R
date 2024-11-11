### cox regression

library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)
library(ggrepel)
library(scales)
library(survival)
library(survminer)




setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

c1c2 <-read_excel("Association/new_pheno_cov/c1c2_eplet_forAsso_withSMMS_1147.xlsx")
head(c1c2)
pheno <- read_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso.xlsx")
head(pheno)


c1c2 %>% select(RecInfo,DonInfo,theme,Total_eps,Total_Abv,Total_Oth,Total_eps_ClassI,Abv_ClassI,Oth_ClassI,Total_eps_ClassII,Abv_ClassII,Oth_ClassII) %>%
  merge(pheno) %>%
  mutate(rej_time_day = ifelse(is.na(rej_time_day),0,rej_time_day))-> c1c2_forAsso
#c1c2_forAsso <- merge(c1c2,pheno)
head(c1c2_forAsso)
dim(c1c2_forAsso)

colnames(c1c2)
c1c2 %>% merge(pheno) %>% filter(theme == "SMMS") %>% #select(Cox_group,DRB,DQB)
  mutate('Risk' = ifelse(DRB ==0 & DQB ==0,"Low", 
                         ifelse(DRB>=7 | DQB >= 9,"High","Intermediate"))) -> c1c2_forCOX_SMMS

c1c2_forCOX_SMMS %>% count(rej_tot,Risk)

coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~Risk,data = c1c2_forCOX_SMMS)
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3)

ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           legend.labs = c("High", "Intermediate","Low"), legend.title = "SM DR-DQ range",) + 
  xlab("Month")


c1c2 %>% merge(pheno) %>% #filter(theme != "SMMS") %>% #select(Cox_group,DRB,DQB)
  select(rej_time_day,theme,rej_tot,Total_eps,Total_eps_ClassI,Total_eps_ClassII,Abv_ClassI,Abv_ClassII,Total_Abv)  %>%
  pivot_longer(cols = 4:9) %>%  #head()
  group_by(theme,name) %>%
  #reframe(q = quantile(value, c(0.25, 0.5, 0.75))) -> c1c2_quantile
  mutate(q=ntile(value,4)) %>% 
  group_by(theme,name,q) %>% #count(q)
  mutate(q_m = max(value)) %>% #head()    
  mutate(q_m = ifelse(q %in% c(2,4),min(value),q_m)) %>% #count(q)
  mutate("g" = ifelse(q==1,"Low",ifelse(q==2,"Intermediate","High"))) -> c1c2_Total_forCox

head(c1c2_Total_forCox)

coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme=="EMS",name == "Total_eps"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "Total_eps") + 
  xlab("Month")


coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme=="EMS",name == "Total_eps_ClassI"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "Total_eps_ClassI") + 
  xlab("Month")


coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme=="EMS",name == "Total_eps_ClassII"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "Total_eps_ClassII") + 
  xlab("Month")

table(c1c2_Total_forCox$name)

coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme=="EMS",name == "Total_Abv"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "Total_Abv") + 
  xlab("Month")


coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme=="EMS",name == "Abv_ClassI"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "Abv_ClassI") + 
  xlab("Month")


coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme=="EMS",name == "Abv_ClassII"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "Abv_ClassII") + 
  xlab("Month")



coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme !="EMS",name == "Total_eps"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "SMMS : Total_eps") + 
  xlab("Month")



coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~g,data = c1c2_Total_forCox %>% filter(theme=="EMS",name == "Total_eps_ClassII"))
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "SMMS : Total_eps_ClassII") + 
  xlab("Month")




table(c1c2_Total_forCox$name)




'''
5.26.(금) 박금보래 교수님 메일
자료 내용을 보고 한가지만 말씀드리면, 
Single molecular mismatch의 경우 기존 논문에서는 DR 0 and DQ 0,  DR 1-6 and/or DQ 1-8,  DR>=7 or DQ>= 9 의 세군으로 분류하고, 
세번째 군을 다시 DQ<15, DQ>=15 로 분류하여 low, intermediate, high risk 세 그룹으로 분석하고  validation 논문까지 나왔습니다.
저희가 weig했을 때 통계적으로  의미를 찾기 어렵다면, 기존연구처럼 아래의 3군으로만 그룹을 나누는 것에 대해서 
한번 고민해 볼 수 있을 것 같습니다.
 
1. low DR 0 and DQ 0  +  DR 1-6 and/or DQ 1-8 (2조건이 모두 low) 
2. intermediate   DR>=7 or DQ>= 9, and  DQ<15
3. high                 DR>=7 or DQ>= 9, and DQ>=15 
'''

'''
myfit <- coxph(Surv(time, event) ~ marker + age + sex + cov1 + cov2, data=data)  # Surv object를 위해 time과 Y값(rej_tot)를 넣고 뒤쪽에 marker(eplet, AMS 등)와 covariate를 넣으면 됨
		summary(myfit)  # 위의 coxph 결과 summary
		summary(myfit)$coefficients  # coeff 추출
		coxsurv <- survfit(myfit)   # survival curve 계산
		ggsurvplot(coxsurv, risk.table=T, fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T, break.time.by=200, conf.int.style="ribbon", conf.int.alpha=0.3)   # survplot 그리는 코드
'''
head(a)
coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = c1c2_forCOX_SMMS)
coxfit <- survfit(coxfit)
coxfit <- ggsurvfit(coxfit)
ggsurvplot(coxfit, risk.table=T, fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T, break.time.by=200, conf.int.style="ribbon", conf.int.alpha=0.3)
head(c1c2_forCOX_SMMS)
table(c1c2_forCOX_SMMS$theme)
coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = c1c2_forCOX_SMMS)
colnames(c1c2_Total_forCox)
coxfit <- survfit(myfit)
summa

ggadjustedcurves(myfit, data= c1c2_forCOX_SMMS, variable="")
