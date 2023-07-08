### new pheno and covariate %20230706
### new sample
library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)
library(ggrepel)
library(scales)



setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")


pheno <- read_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso.xlsx")
pheno <- read_xlsx("00.pheno/20221109_KOTRY_selectFeature_forAuteRejectionAsso_v2.xlsx")
pheno %>% mutate(rej_time_day = date_rej_tot - date_kt) ->pheno
plot(pheno$rej_time_day,pheno$rej_tot)
table(pheno$rej_tot)
head(pheno)
dim(pheno)
ref <- read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt",header = T)

c1 <-read_xlsx("Association/c1_forWAS.xlsx")
c2 <-read_xlsx("Association/c2_forWAS.xlsx")


c1c2 <-read_xlsx("Association/new_pheno_cov/c1c2_eplet_forAsso_withSMMS_1147.xlsx")

colnames(c1c2)
dim(c1c2)
#c1c2 %>% select(1:28)
c1_c2_eplet_forAsso <-read_xlsx("Association/new_pheno_cov/c1c2_eplet_forAsso_withSMMS_1147.xlsx") %>% merge(pheno)
c1_c2_eplet_forAsso %>% filter(theme == "EMS") %>% count(rej_tot)

plot(c1_c2_eplet_forAsso$rej_time_day,c1_c2_eplet_forAsso$rej_tot)


head(c1_c2_eplet_forAsso)
#colnames(pheno)

colnames(c1_c2_eplet_forAsso)
dim(c1_c2_eplet_forAsso)
result <- NULL
for(i in 3:27){
  if (colnames(c1_c2_eplet_forAsso)[i] != "theme") {
    c1_c2_eplet_forAsso[,i] <- scale(as.numeric(c1_c2_eplet_forAsso[,i]))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme == "EMS"), family="binomial")
    result <- rbind(result, c("EMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
    temp <- glm(paste("rej_tot ~ ",colnames(c1_c2_eplet_forAsso)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c1_c2_eplet_forAsso %>% filter(theme != "EMS"), family="binomial")
    result <- rbind(result, c("SMMS",colnames(c1_c2_eplet_forAsso)[i],as.vector(summary(temp)$coefficients[2,])))
  }  
}
head(result)
as.vector(summary(temp)$coefficients[2,])
summary(temp)

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
filter(!(theme == "SMMS" & !(ID %in% c("Total_eps","Total_Abv","Total_Oth","Total_eps_ClassII","Abv_ClassII","Oth_ClassII","AbDQB","AbDRB","DQB","DRB","otDQB","otDRB")))) -> result_forplot


#writexl::write_xlsx(result_forplot,"Association/new_pheno_cov/AuteRejection_newPheno20221109_c1c2_eplet_assoResult_withouthHLAms_scale_withSMMS_onlyEpletCount.xlsx")

df <- read_xlsx("Association/new_pheno_cov/AuteRejection_newPheno20221109_c1c2_eplet_assoResult_withouthHLAms_scale_withSMMS_onlyEpletCount.xlsx")
head(df)
table(df$ID)
df %>% mutate(ID = ifelse(grepl("Total_eps",ID),"ALL",
                          ifelse(grepl("Abv",ID),"Abv",
                                 ifelse(grepl("Oth",ID),"Oth",ID)))) -> df




table(c1_c2_eplet_forAsso$rej_tot)/2
sum(table(c1_c2_eplet_forAsso$rej_tot)/2)

head(df)
#ggplot(aes(x=ID,y=-log10P,color=factor(class,levels = c("CLASS I + II","CLASS I","CLASS II")),shape=theme,size=Estimate)) + 
ggplot(df,aes(x=ID,y=-log10P,shape=theme,colour=Estimate)) + 
  geom_point(size = 2) +
  scale_color_gradient2(low="blue",mid ="green",high="red",midpoint =0.1) +
  labs(colour='Estimate',
       shape = "")  +
  theme(axis.text.x = element_text(size = 10),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        #legend.title = element_blank(),
        legend.text = element_text(size = 13)) + 
  theme(strip.text.x = element_text(size = 13,face = "bold")) + 
  facet_grid(~factor(class,levels = c("CLASS I + II","CLASS I","CLASS II","CLASS II (Gene)")),scales = "free_x",space = "free_x") + 
  geom_text_repel(#data=subset(df,-log10P > 1.25),
                  data=subset(df,P < 0.05),
                  aes(label=ID),
                  show.legend = F,
                  size =4,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))




####### Amino Acid #####################################################################################################################################
head(c1)
head(pheno)

pheno %>% merge(c1) -> c1_data


head(c1_data)
table(c1_data$rej_tot)
colnames(c1_data)[30:ncol(c1_data)]
colnames(c1_data)[45:ncol(c1_data)]

colnames(c1_data)[48:ncol(c1_data)] <- paste0("AA_",colnames(c1_data)[48:ncol(c1_data)])

#c1_data_scale <- c1_data

for (i in 48:ncol(c1_data)) {
  c1_data[,i] <- as.numeric(c1_data[,i])
}


result <- NULL
#result <- NULL
for(i in 48:ncol(c1_data)){
  if (length(table(c1_data[,i])) == 2 ) {
    temp <- glm(paste("rej_tot ~ ",colnames(c1_data)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c1_data, family="binomial")
    result <- rbind(result, c(colnames(c1_data)[i],as.vector(summary(temp)$coefficients[2,])))
  }
}


head(result)
colnames(result) <- c("ID","Estimate","SE","Z","P")

c1_result <- result %>% as.data.frame()
c1_result

#write.table(c1_result, "Association/new_pheno_cov/c1_aplet_AAWAS_newpheno.txt", col.names=T, row.names=F, sep="\t", quote=F)


## 253Q --Abv
## 11AV --other
#c1_result <- read_table("Association/c1_eplet_asso.txt")
c1_result <- read_table("Association/new_pheno_cov/c1_aplet_AAWAS_newpheno.txt")
head(c1_result)
head(c1_data)
colnames(c1_result)
c1_result$ID
c1_result$ID[65:75]
str_replace(c1_result$ID,"AA_","AbV_")
c1_result$ID[1:67] <- str_replace(c1_result$ID[1:67],"AA_","AbV_")
c1_result$ID[68:nrow(c1_result)] <- str_replace(c1_result$ID[68:nrow(c1_result)],"AA_","Oth_")

head(c1_result)

head(c1_result)
c1_result %>% as.data.frame() %>% mutate('log10P' = log10(as.numeric(P))) %>% filter(!ID %in% c("AbVer_Eplets","Other_Eplets") ) %>%
  mutate(type= ifelse(grepl("AbV_",ID),"Amino Acid (AbV)",ifelse(grepl("Oth_",ID),"Amino Acid (Other)","Eplet_Count"))) %>% 
  mutate(Estimate = as.numeric(Estimate))-> c1_result
#mutate_all(funs(str_replace(.,"AA_",""))) %>% head()
head(c1_result)


ggplot(c1_result,aes(x=fct_inorder(ID),y = -log10P,color = type)) +
  geom_point() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
  geom_text_repel(data=subset(c1_result,P < 0.05),
                  aes(label=ID),
                  show.legend = F,
                  size =5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))

c1_result %>% filter(P < 0.05) -> a
a
head(c1_result) -> a


###c2#######
pheno %>% merge(c2) -> c2_data

head(c2_data)
table(c2_data$rej_tot)
colnames(c2_data)[43:ncol(c2_data)]
colnames(c2_data)[60:ncol(c2_data)]

colnames(c2_data) %>% as_data_frame() %>% count(value) %>% filter(n==2)
colnames(c2_data)[i]

#colnames(data2)[36:ncol(data2)] <- paste0("AA_",colnames(data2)[36:ncol(data2)])

colnames(c2_data)
colnames(c2_data)[75:85]

for (i in 81:ncol(c2_data)) {
  c2_data[,i] <- as.numeric(c2_data[,i])
  
}


head(c2_data)
table(c2_data[,i])
length(table(c2_data[,i]))

result <- NULL

for(i in 81:ncol(c2_data)){
  if (length(table(c2_data[,i])) == 2 ) {
    temp <- glm(paste("rej_tot ~ ",colnames(c2_data)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group",sep=""), data=c2_data, family="binomial")
    result <- rbind(result, c(colnames(c2_data)[i],as.vector(summary(temp)$coefficients[2,])))
  }
}

colnames(result) <- c("ID","Estimate","SE","Z","P")

result %>% as.data.frame()


c2_result <- result %>% as.data.frame()

#write.table(c2_result, "Association/new_pheno_cov/c2_aplet_AAWAS_newpheno.txt", col.names=T, row.names=F, sep="\t", quote=F)


c2_result <- read_table("Association/new_pheno_cov/c2_aplet_AAWAS_newpheno.txt")

head(c2_result)
grep("DR",c2_result$ID)
grepl("DR",c2_result$ID)
c2_result$ID
c2_result %>% as.data.frame() %>% mutate('log10P' = log10(as.numeric(P))) %>% #max(log10P)
  mutate('type' = ifelse(grepl("DR",ID),"DR",ifelse(grepl("DQ",ID),"DQ",ifelse(grepl("DP",ID),"DP","Total")))) -> c2_p


c2_result %>% as.data.frame() %>% mutate('log10P' = log10(as.numeric(P))) %>% #max(log10P)
  mutate('type' = ifelse(grepl("DR",ID),"DRB1",ifelse(grepl("DQB",ID),"DQB1",ifelse(grepl("DP",ID),"DPB1",ifelse(grepl("QA",ID),"DQA1",ifelse(grepl("PA",ID),"DPA1","check")))))) %>% 
  ggplot(aes(x=ID,y = -log10P,color = type)) +
  geom_point() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

c2_result$ID

c2_result %>% filter(ID %in% c("All_Eplet","All_AbV","Nr_allDRB","Nr_AbDRB","Nr_otDRB","Nr_allDQB",
                               "Nr_AbDQB","Nr_otDQB","Nr_allDQA","Nr_AbDQA","Nr_otDQA","Nr_allDPB","Nr_AbDPB","Nr_otDPB",
                               "Nr_allDPA","Nr_AbDPA","Nr_otDPA")) %>% #head()
  mutate('log10P' = log10(as.numeric(P))) %>% #max(log10P)
  mutate('type' = ifelse(grepl("DR",ID),"DR",ifelse(grepl("DQ",ID),"DQ",ifelse(grepl("DP",ID),"DP","Total")))) %>% 
  ggplot(aes(x=fct_inorder(ID),y=-log10P,color=type)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        legend.title = element_blank())
#scale_x_discrete(breaks=arrange(., ID)$id, labels=arrange(., ID)$label))


c2_result$ID


c2_result %>% filter(!ID %in% c("All_AbV","All_Eplet","Nr_allDRB","Nr_AbDRB","Nr_otDRB","Nr_allDQB",
                                "Nr_AbDQB","Nr_otDQB","Nr_allDQA","Nr_AbDQA","Nr_otDQA","Nr_allDPB","Nr_AbDPB","Nr_otDPB",
                                "Nr_allDPA","Nr_AbDPA","Nr_otDPA")) %>% #head()
  mutate('log10P' = log10(as.numeric(P))) %>% #head()
  #filter(log10P > -1.5) %>%
  mutate('Gene' = ifelse(grepl("DR",ID),"DRB1",ifelse(grepl("DQB",ID),"DQB1",ifelse(grepl("DP",ID),"DPB1",ifelse(grepl("QA",ID),"DQA1",ifelse(grepl("PA",ID),"DPA1","check")))))) %>% 
  mutate('type' = ifelse(grepl("Ab",ID),"AbV",ifelse(grepl("Ot",ID) | grepl("ot",ID) ,"Other","E"))) -> c2_result_v2


head(c2_result_v2)

ggplot(c2_result_v2,aes(x=fct_inorder(ID),y=-log10P,color=fct_inorder(Gene),shape=type)) +
  geom_point(size = 2) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom") + 
  geom_text_repel(data=subset(c2_result_v2,P < 0.05),
                  aes(label=ID),
                  show.legend = F,
                  size =5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))

c2_result_v2 %>% filter(P < 0.05) -> a

'''
5.26.(금) 박금보래 교수님 메일
자료 내용을 보고 한가지만 말씀드리면, 
Single molecular mismatch의 경우 기존 논문에서는 DR 0 and DQ 0,  DR 1-6 and/or DQ 1-8,  DR>=7 or DQ>= 9 의 세군으로 분류하고, 
세번째 군을 다시 DQ<15, DQ>=15 로 분류하여 low, intermediate, high risk 세 그룹으로 분석하고  validation 논문까지 나왔습니다.
저희가 weighted로 분석했을 때 통계적으로  의미를 찾기 어렵다면, 기존연구처럼 아래의 3군으로만 그룹을 나누는 것에 대해서 
한번 고민해 볼 수 있을 것 같습니다.
 
1. low DR 0 and DQ 0  +  DR 1-6 and/or DQ 1-8 (2조건이 모두 low) 
2. intermediate   DR>=7 or DQ>= 9, and  DQ<15
3. high                 DR>=7 or DQ>= 9, and DQ>=15 
'''


###########

### cox
#coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% filter(theme == "EMS"))
tmp <- coxph(Surv(rej_time_day/30,rej_tot)~DRB+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = c1_c2_eplet_forAsso %>% filter(theme == "EMS"))
head(c1_c2_eplet_forAsso)
tmp
result <- NULL
head(c1_c2_eplet_forAsso)
#c1_c2_eplet_forAsso %>% filter(theme == "EMS")
c1_c2_eplet_forAsso %>% filter(theme == "EMS") %>% as.data.frame()-> a
colnames(c1_c2_eplet_forAsso)
for(i in 3:27){
  if (colnames(a)[i] != "theme") {
    a[,i] <- scale(as.numeric(a[,i]))
    print(colnames(a)[i])
    my_formula <- as.formula(paste0("Surv(rej_time_day, rej_tot)~",colnames(a)[i],"+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group"))
    print(my_formula)
    temp <- coxph(my_formula,data=a)
    result <- rbind(result,c("EMS",colnames(a)[i],as.vector(summary(temp)$coefficients[1,])))
  }
}


### cox
#coxfit <- coxph(Surv(rej_time_day/30,rej_tot)~Risk+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = tmp %>% filter(theme == "EMS"))
tmp <- coxph(Surv(rej_time_day/30,rej_tot)~DRB+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group,data = c1_c2_eplet_forAsso %>% filter(theme == "EMS"))
head(c1_c2_eplet_forAsso)
tmp
result <- NULL
head(c1_c2_eplet_forAsso)
#c1_c2_eplet_forAsso %>% filter(theme == "EMS")
c1_c2_eplet_forAsso %>% filter(theme == "EMS") %>% as.data.frame()-> a
head(c1_c2_eplet_forAsso)
colnames(c1_c2_eplet_forAsso)
for(i in 3:27){
  if (colnames(a)[i] != "theme") {
    a[,i] <- scale(as.numeric(a[,i]))
    print(colnames(a)[i])
    my_formula <- as.formula(paste0("Surv(rej_time_day, rej_tot)~",colnames(a)[i],"+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf+TACDR+CYCDR+desen+group"))
    print(my_formula)
    temp <- coxph(my_formula,data=a)
    result <- rbind(result,c("EMS",colnames(a)[i],as.vector(summary(temp)$coefficients[1,])))
  }
}
result
head(a)

a %>% count(rej_tot)
plot(a$rej_time_day,a$rej_tot)
ggforest(temp, data=a)

head(result)

result <- result %>% as.data.frame()
colnames(result) <- c("Method","Count","coef","exp(coef)","se(coef)","Z","P")
head(ref)
table(result$Method)
#writexl::write_xlsx(result,"./Association/new_pheno_cov/CoxResult.eplet.count_withSMMS_notscale.xlsx")
