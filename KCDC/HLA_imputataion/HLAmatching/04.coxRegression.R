### cox regression

library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)
library(ggrepel)
library(scales)
library(survival)
library(survminer)
library(cowplot)



setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

c1c2 <-read_excel("Association/new_pheno_cov/c1c2_eplet_forAsso_withSMMS_1147.xlsx")
head(c1c2)
pheno <- read_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso.xlsx")
pheno <- read_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso_v2.xlsx")
head(pheno)

pheno %>% mutate(rej_time_day = date_rej_tot - date_kt) ->pheno
#pheno %>% mutate(rej_time_day = date_rej_tot - date_kt) %>% select(rej_time_day)

c1c2 %>% select(RecInfo,DonInfo,theme,Total_eps,Total_Abv,Total_Oth,Total_eps_ClassI,Abv_ClassI,Oth_ClassI,Total_eps_ClassII,Abv_ClassII,Oth_ClassII) %>%
  merge(pheno) %>%
  mutate(rej_time_day = ifelse(is.na(rej_time_day),0,rej_time_day))-> c1c2_forAsso

plot(c1c2_forAsso$rej_time_day,c1c2_forAsso$rej_tot)
#c1c2_forAsso <- merge(c1c2,pheno)
head(c1c2_forAsso)
head(c1c2_forAsso$rej_time_day)
table(is.na(c1c2_forAsso$rej_time_day))
dim(c1c2_forAsso)

colnames(c1c2)
c1c2 %>% merge(pheno) %>% filter(theme == "SMMS") %>% #select(Cox_group,DRB,DQB)
  mutate('Risk' = ifelse(DRB ==0 & DQB ==0,0, 
                     ifelse(DRB>=7 | DQB >= 9,2,1))) -> c1c2_forCOX_SMMS

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
  #mutate("g" = ifelse(q==1,"Low",ifelse(q==4,"High","Intermediate"))) -> c1c2_Total_forCox
  mutate("g" = ifelse(q==1,0,ifelse(q==4,2,1))) -> c1c2_Total_forCox

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


''' test



'''
5.26.(금) 박금보래 교수님 메일
자료 내용을 보고 한가지만 말씀드리면, 
Single molecular mismatch의 경우 기존 논문에서는 DR 0 and DQ 0,  DR 1-6 and/or DQ 1-8,  DR>=7 or DQ>= 9 의 세군으로 분류하고, 
세번째 군을 다시 DQ<15, DQ>=15 로 분류하여 low, intermediate, high risk 세 그룹으로 분석하고  validation 논문까지 나왔습니다.
저희가 weighted로 분석했을 때 통계적으로  의미를 찾기 어렵다면, 기존연구처럼 아래의 3군으로만 그룹을 나누는 것에 대해서 
한번 고민해 볼 수 있을 것 같습니다.
''' 
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

head(c1c2_forCOX_SMMS)
table(c1c2_forCOX_SMMS$theme)
table(c1c2_forCOX_SMMS$Risk)
coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = c1c2_forCOX_SMMS)
summary(coxfit)

ggadjustedcurves(coxfit, data=c1c2_forCOX_SMMS, variable="Risk",xlab="Month",legend.title = "SMMS DR-DQ Range")

head(c1c2_forAsso)
#table(pheno$group)
c1c2_forAsso %>% 
  group_by(theme) %>% 
  mutate(q=ntile(Total_eps,4)) %>% #count(q)
  #mutate("Risk" = ifelse(q==1,"Low",ifelse(q==4,"High","Intermediate"))) -> c1c2_forAsso_cox
  mutate("Risk" = ifelse(q==1,0,ifelse(q==4,2,1))) -> c1c2_forAsso_cox

colnames(c1c2_forAsso_cox)
table(c1c2_forAsso_cox$Risk)
table(c1c2_forAsso_cox$group)
c1c2_forAsso_cox %>% filter(theme == "EMS") %>% as.data.frame() -> a
head(a)
head(c1c2_forAsso_cox %>% filter(theme == "EMS"))
head(c1c2_forAsso_cox$rej_time_day)
str(c1c2_forAsso_cox$rej_time_day)

result <- NULL

head(c1c2_forAsso_cox)
head(c1c2_forAsso)

c1c2_forAsso %>%
  plot(rej_time_day,rej_tot)
plot(c1c2_forAsso$rej_time_day,c1c2_forAsso$rej_tot)
#20201231 rejet tot == 0

summary(coxfit)$coefficients[1,]
coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = c1c2_forAsso_cox %>% filter(theme == "EMS") %>% as.data.frame())
ggadjustedcurves(coxfit, data=c1c2_forAsso_cox %>% filter(theme == "EMS") %>% as.data.frame(), variable="Risk",xlab="Month",legend.title = "Total eps (EMS)")  -> a

result <- rbind(result,c("EMS","Total_eps","basic",as.vector(summary(coxfit)$coefficients[1,])))
result

coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = c1c2_forAsso_cox %>% filter(theme != "EMS") %>% as.data.frame())
ggadjustedcurves(coxfit, data=c1c2_forAsso_cox %>% filter(theme != "EMS") %>% as.data.frame(), variable="Risk",xlab="Month",legend.title = "Total eps (SMMS)")  -> b

result <- rbind(result,c("SMMS","Total_eps","basic",as.vector(summary(coxfit)$coefficients[1,])))
result
plot_grid(a,b,labels = c("EMS","SMMS"))

head(c1c2_forAsso)
c1c2_forAsso %>% as_tibble() %>%
  mutate(date_rej_tot = ifelse(is.na(date_rej_tot),"20201231",date_rej_tot)) %>% select(date_rej_tot)
#20201231 rejet tot == 0
date("2020-12-31")
date("20201231")
colnames(c1c2_forAsso)
c1c2_forAsso %>% 
  group_by(theme) %>% 
  mutate(q=ntile(Total_eps,4)) %>% #count(q)
  mutate("Risk" = ifelse(q==1,0,ifelse(q==4,2,1))) %>%
  as.data.frame() -> tmp

coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% filter(theme == "EMS"))
ggadjustedcurves(coxfit, data=tmp %>% filter(theme == "EMS"),  variable="Risk",xlab="Month",legend.title = "Total eps (EMS)")  -> a

result <- rbind(result,c("EMS","Total_eps","basic",as.vector(summary(coxfit)$coefficients[1,])))
result

colnames(c1c2_forAsso)
c1c2_forAsso %>% 
  group_by(theme) %>% 
  mutate(q=ntile(Total_eps_ClassI,4)) %>% #count(q)
  mutate("Risk" = ifelse(q==1,0,ifelse(q==4,2,1))) %>%
  as.data.frame() -> tmp

coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% filter(theme == "EMS"))
ggadjustedcurves(coxfit, data=tmp %>% filter(theme == "EMS"),  variable="Risk",xlab="Month",legend.title = "Total eps ClassI (EMS)")  -> b

result <- rbind(result,c("EMS","Total_eps_ClassI","basic",as.vector(summary(coxfit)$coefficients[1,])))
result
colnames(c1c2_forAsso)
c1c2_forAsso %>% 
  group_by(theme) %>% 
  mutate(q=ntile(Total_eps_ClassII,4)) %>% #count(q)
  mutate("Risk" = ifelse(q==1,0,ifelse(q==4,2,1))) %>%
  as.data.frame() -> tmp

coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% filter(theme == "EMS"))
ggadjustedcurves(coxfit, data=tmp %>% filter(theme == "EMS"),  variable="Risk",xlab="Month",legend.title = "Total eps ClassII (EMS)")  -> c
result <- rbind(result,c("EMS","Total_eps_ClassII","basic",as.vector(summary(coxfit)$coefficients[1,])))
result
plot_grid(a,b,c, labels = c("Total EPS","Total EPS ClassI","Total EPS ClassII"),ncol = 3)





colnames(c1c2_forAsso)
c1c2_forAsso %>% 
  group_by(theme) %>% 
  mutate(q=ntile(Total_eps,4)) %>% #count(q)
  mutate("Risk" = ifelse(q==1,0,ifelse(q==4,2,1))) %>%
  as.data.frame() -> tmp

coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% filter(theme != "EMS"))
ggadjustedcurves(coxfit, data=tmp %>% filter(theme != "EMS"),  variable="Risk",xlab="Month",legend.title = "Total eps (SMMS)")  -> a
result <- rbind(result,c("SMMS","Total_eps","basic",as.vector(summary(coxfit)$coefficients[1,])))
result
colnames(c1c2_forAsso)
c1c2_forAsso %>% 
  group_by(theme) %>% 
  mutate(q=ntile(Total_eps_ClassII,4)) %>% #count(q)
  mutate("Risk" = ifelse(q==1,0,ifelse(q==4,2,1))) %>%
  as.data.frame() -> tmp

coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% filter(theme != "EMS"))
ggadjustedcurves(coxfit, data=tmp %>% filter(theme != "EMS"),  variable="Risk",xlab="Month",legend.title = "Total eps ClassII (SMMS)")  -> b
result <- rbind(result,c("SMMS","Total_eps_ClaaII","basic",as.vector(summary(coxfit)$coefficients[1,])))
result
c1c2 %>% merge(pheno) %>% filter(theme == "SMMS") %>% 
  mutate('Risk' = ifelse(DRB ==0 & DQB ==0,0,
                  ifelse((DRB >= 1 & DRB <= 6) |(DQB >= 1 & DRB <= 8),0,
                  ifelse((DRB>=7 | DQB>=9) & DQB >= 15),2,1)))) -> tmp
  #mutate('Risk' = ifelse(DRB ==0 & DQB ==0,0, 
                         #ifelse(DRB>=7 | DQB >= 9,2,1))) -> tmp
tmp %>% count(Risk)

'''
1. low DR 0 and DQ 0  +  DR 1-6 and/or DQ 1-8 (2조건이 모두 low) 
2. intermediate   DR>=7 or DQ>= 9, and  DQ<15
3. high                 DR>=7 or DQ>= 9, and DQ>=15 
'''

#coxfit <- coxph(Surv(month(rej_time_day),rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp)
coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp)
ggadjustedcurves(coxfit, data=tmp,  variable="Risk",xlab="Month",legend.title = "DR-DQ (SMMS)")  -> c
result <- rbind(result,c("SMMS","DR-DQ","byPark",as.vector(summary(coxfit)$coefficients[1,])))
result
plot_grid(a,b,c, labels = c("A : Total EPS","B : Total EPS ClassII","C :DR-DQ"),ncol = 3)
c

tmp %>% head()
hist(as.numeric(tmp$rej_time_day/30))

tmp %>% ggplot(aes(x=rej_time_day/30)) + 
  geom_histogram() + 
  facet_grid(~rej_tot)

tmp$rej_time_day/30
month(tmp$rej_time_day)

result
summary(coxfit)


tmp %>% count(rej_tot,rej_time_day)


c1c2_forAsso %>% 
  ggplot(aes(x=rej_time_day)) + 
  geom_histogram() + 
  facet_grid(~rej_tot)


coxfit <- coxph(Surv(rej_time_day,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp)
ggadjustedcurves(coxfit, data=tmp,  variable="Risk",xlab="day",legend.title = "DR-DQ (SMMS)")

coxfit <- (Surv(rej_time_day,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp)
ggadjustedcurves(coxfit, data=tmp,  variable="Risk",xlab="day",legend.title = "DR-DQ (SMMS)")

coxfit <- survfit(Surv(rej_time_day/30,rej_tot)~Risk,data = tmp)
ggsurvplot(coxfit, risk.table=T,fun='pct', pval=T,pval.coord=c(0.2,0.2), conf.int=T,conf.int.style="ribbon", conf.int.alpha=0.3,
           #           legend.labs = c("High", "Intermediate","Low"), legend.title = "Total_eps") +
           legend.title = "DR-DQ (SMMS)") + 
  xlab("Month")
