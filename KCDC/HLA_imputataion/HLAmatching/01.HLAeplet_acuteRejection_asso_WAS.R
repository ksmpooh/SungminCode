## HLA eplet aa GWAS
library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)
library(ggrepel)

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230508_epletv2/")


c1_header <- read_xlsb("ABC_Eplet_Matching_4.0_g1.xlsb",sheet = 4,skip = 2)

c1_header[1,]
c1_1 <- read_xlsb("ABC_Eplet_Matching_4.0_g1.xlsb",sheet = 4,skip = 4)
c1_2 <- read_xlsb("ABC_Eplet_Matching_4.0_g2.xlsb",sheet = 4,skip = 4)
colnames(c1_1)
#df %>% select(grep("column",colnames(df)),X,X.1,X.2,X.3,X.4,X.5,X.6) 
c1_1[,c(1:20)] %>% #tail()
  filter(X != "x") %>% filter(!is.na(X)) %>% #tail()  
  filter(X != "") %>%  head()


#colnames(c1_1) <- c1_header[1,]
#colnames(c1_2) <- c1_header[1,]
head(c1_1)
head(c1_2)

c1_1 <- c1_1 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")
c1_2 <- c1_2 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")

c1_header[1,]
c1_header[1,] %>% t() %>% as.data.frame() %>% #head()
  rename("col" = 1) %>% count(col)
colnames(c1_1)
#c1_1 %>% rbind(c1_2) %>% 

'
c1_2[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-") %>% tail()

c1_1 <- c1_1[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")

c1_1 <- c1_1[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")
c1_2 <- c1_2[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")
tail(c1_1)
head(c1_1)
'
#c1_header[1,] %>% write_xlsx("c1_header.xlsx")
# 이거 후에 수작업


c1_hh <- read_excel("c1_header.xlsx",sheet = 2)
c1_hh
c1_hh %>% t() %>% as.data.frame() -> c1_header_index 

head(c1_header_index)
c1 <- rbind(c1_1,c1_2)
head(c1)
seq(1:ncol(c1))
colnames(c1) <- seq(1:ncol(c1))

head(c1_header_index)
c1 %>% select(c1_header_index$V2) -> c1
colnames(c1) <- c1_header_index$V1
head(c1)
head(c1[,20:26])
colnames(c1)[20:26]
c1 %>% mutate_all(~ replace(.,is.na(.),0)) %>%
  mutate(across(colnames(c1)[21:ncol(c1)], ~ replace(.,.!=0 , 1))) -> c1_forWAS




c2_1 <- read_xlsb("DRDQDP_Eplet_Matching_3.1_g1.xlsb",sheet = 4,skip = 4)
c2_2 <- read_xlsb("DRDQDP_Eplet_Matching_3.1_g2.xlsb",sheet = 4,skip = 4)
colnames(c2_2)

head(c2_1)
head(c2_2)


c2_1 <- c2_1 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")
c2_2 <- c2_2 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")

head(c2_1)
head(c2_2)
c2 <- rbind(c2_1,c2_2)
head(c2)

'
c2_header <- read_xlsb("DRDQDP_Eplet_Matching_3.1.xlsb",sheet = 4,skip = 1)
#c2_header[1,] %>% t()

c2_header[1,] %>% write_xlsx("c2_header_extra.xlsx")
# 이거 후에 수작업
'

c2_hh <- read_excel("c2_header.xlsx",sheet = 4)
c2_hh
#c2_hh %>% t() %>% as.data.frame() -> c2_header_index 

head(c2_hh)
c2 <- rbind(c2_1,c2_2)
head(c2)
seq(1:ncol(c2))
colnames(c2) <- seq(1:ncol(c2))

head(c2)
c2 %>% select(c2_hh$V4) -> c2
colnames(c2) <- c2_hh$V5
head(c2)
head(c2[,52:80])
colnames(c2)[52:80]
c2 %>% mutate_all(~ replace(.,is.na(.),0)) %>%
  mutate(across(colnames(c2)[55:ncol(c2)], ~ replace(.,.!=0 , 1))) -> c2_forWAS

colnames(c2) %>% as_data_frame() %>% count(value) %>% filter(n != 1)

head(c2_forWAS)
head(c1_forWAS)


#write_xlsx(c1_forWAS,"~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/Association/c1_forWAS.xlsx")
#write_xlsx(c2_forWAS,"~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/Association/c2_forWAS.xlsx")


############

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/Association/")

pheno <- read_table("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/00.pheno/Rejection_phenotype_coded_KID_20230323.txt")
ref <- read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt",header = T)

c1 <- read_xlsx("c1_forWAS.xlsx")
c2 <- read_xlsx("c2_forWAS.xlsx")
head(pheno)
head(ref)
head(c1)
head(c2)



#### eplet histogram
c1 %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>%
  pivot_longer(cols = 3:5,names_to = "type",values_to = "score") %>% #head()
  ggplot(aes(x=score,fill=type))+
  geom_histogram() +
  facet_grid(~fct_inorder(type)) +
  theme(legend.position = 'none',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 13,face = "bold"))

colnames(c2)
c2 %>% select(RecInfo,DonInfo,28:44) %>% colnames()
c2 %>% select(RecInfo,DonInfo,28:44) %>% #head()#dim()
  pivot_longer(cols = 3:4,names_to = "type",values_to = "score") %>% #head()
  ggplot(aes(x=score,fill=type))+
  geom_histogram() +
  facet_grid(~type) +
  theme(legend.position = 'none',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 13,face = "bold"))

c2 %>% select(RecInfo,DonInfo,28:44) %>%# head()#dim()
  pivot_longer(cols = c(5,8,11,14,17),names_to = "type",values_to = "score") %>%
  ggplot(aes(x=score,fill=type))+
  geom_histogram() +
  facet_grid(~type) +
  theme(legend.position = 'none',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 13,face = "bold"))


###################################


pheno %>% merge(ref,by.x = "KCHIP_ID",by.y="KBA_ID") %>% 
  merge(c1,by.x="bCODE",by.y="RecInfo") -> c1_data

head(c1_data)
table(c1_data$rej_tot)
colnames(c1_data)[30:ncol(c1_data)]

colnames(c1_data)[36:ncol(c1_data)] <- paste0("AA_",colnames(c1_data)[36:ncol(c1_data)])

for (i in 36:ncol(c1_data)) {
  c1_data[,i] <- as.numeric(c1_data[,i])
}

result <- NULL
for(i in 31:33){
  temp <- glm(paste("rej_tot ~ ",colnames(c1_data)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c1_data, family="binomial")
  result <- rbind(result, c(colnames(c1_data)[i],as.vector(summary(temp)$coefficients[2,])))
}


#result <- NULL
for(i in 36:ncol(c1_data)){
  if (length(table(c1_data[,i])) == 2 ) {
    temp <- glm(paste("rej_tot ~ ",colnames(c1_data)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c1_data, family="binomial")
    result <- rbind(result, c(colnames(c1_data)[i],as.vector(summary(temp)$coefficients[2,])))
  }
}








head(result)
colnames(result) <- c("ID","Estimate","SE","Z","P")
c1_result <- result %>% as.data.frame()

#write.table(c1_result, "c1_eplet_asso.txt", col.names=T, row.names=F, sep="\t", quote=F)
#write.table(c1_result, "c1_eplet_asso_withouthHLAms.txt", col.names=T, row.names=F, sep="\t", quote=F)


## 253Q --Abv
## 11AV --other
#c1_result <- read_table("c1_eplet_asso.txt")
c1_result <- read_table("c1_eplet_asso_withouthHLAms.txt")
head(c1_result)
head(c1_data)
colnames(c1_result)
c1_result$ID
c1_result$ID[71:75]
str_replace(c1_result$ID,"AA_","AbV_")
c1_result$ID[6:73] <- str_replace(c1_result$ID[6:73],"AA_","AbV_")
c1_result$ID[74:nrow(c1_result)] <- str_replace(c1_result$ID[74:nrow(c1_result)],"AA_","Oth_")

head(c1_result)

c1_result %>% as.data.frame() %>% mutate('-log10P' = -log10(as.numeric(P))) %>% 
  filter(!ID %in% c("AbVer_Eplets","Other_Eplets") ) %>% #head()
  mutate(type= ifelse(grepl("AA",ID),"Amino Acid","Eplet_Count")) %>% # head()
  #mutate_all(funs(str_replace(.,"AA_",""))) %>% head()
  ggplot(aes(x=fct_inorder(ID),y = -log10P,color = type)) +
  geom_point() + 
  geom_text_repel(data=subset(fct_inorder(ID),-log10P > 2.0),
                  aes(label=ID),
                  size =5, 
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom") 


head(c1_result)
c1_result %>% as.data.frame() %>% mutate('log10P' = log10(as.numeric(P))) %>% filter(!ID %in% c("AbVer_Eplets","Other_Eplets") ) %>%
  mutate(type= ifelse(grepl("AbV_",ID),"Amino Acid (AbV)",ifelse(grepl("Oth_",ID),"Amino Acid (Other)","Eplet_Count"))) -> c1_result
  #mutate_all(funs(str_replace(.,"AA_",""))) %>% head()

ggplot(c1_result,aes(x=fct_inorder(ID),y = -log10P,color = type)) +
  geom_point() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom") +
  geom_text_repel(data=subset(c1_result,-log10P > 1.5),
                  aes(label=ID),
                  show.legend = F,
                  size =5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))


c1_result %>% filter(log10P < -1.75) -> a
head(c1_result) -> a
###c2

c2 <- read_xlsx("c2_forWAS.xlsx")

pheno %>% merge(ref,by.x = "KCHIP_ID",by.y="KBA_ID") %>% 
  merge(c2,by.x="bCODE",by.y="RecInfo") %>%
  select(!grep("Eps_",colnames(.))) -> c2_data
  

head(c2_data)
table(c2_data$rej_tot)
colnames(c2_data)[43:ncol(c2_data)]
colnames(c2_data)[60:ncol(c2_data)]

colnames(c2_data) %>% as_data_frame() %>% count(value) %>% filter(n==2)

colnames(c2_data)[i]

#colnames(data2)[36:ncol(data2)] <- paste0("AA_",colnames(data2)[36:ncol(data2)])

str(c2_data)
str(c1_data)

for (i in 60:ncol(c2_data)) {
  c2_data[,i] <- as.numeric(c2_data[,i])
}
head(c2_data)
table(c2_data[,i])
length(table(c2_data[,i]))

head(c2_data[43:45])
head(c2_data[55:60])
head(c2_data[59:65])

result <- NULL
for(i in 43:59){
  #table(c2_data[,i])
  temp <- glm(paste("rej_tot ~ ",colnames(c2_data)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_data, family="binomial")
  result <- rbind(result, c(colnames(c2_data)[i],as.vector(summary(temp)$coefficients[2,])))
}

#result <- NULL
for(i in 60:ncol(c2_data)){
  if (length(table(c2_data[,i])) == 2 ) {
    temp <- glm(paste("rej_tot ~ ",colnames(c2_data)[i], "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_data, family="binomial")
    result <- rbind(result, c(colnames(c2_data)[i],as.vector(summary(temp)$coefficients[2,])))
  }
}



head(result)
colnames(result) <- c("ID","Estimate","SE","Z","P")

result %>% as.data.frame()


c2_result <- result %>% as.data.frame()


c2_result %>% mutate(ID = ifelse(ID=="All_ABV_ClassII","All_AbV",ID)) %>%
  mutate(ID = ifelse(ID=="Total_eps","All_Eplet",ID)) ->  c2_result
  

#write.table(c2_result, "c2_eplet_asso.txt", col.names=T, row.names=F, sep="\t", quote=F)
#write.table(c2_result, "c2_eplet_asso_withouthHLAms.txt", col.names=T, row.names=F, sep="\t", quote=F)

#c2_result <- read_table("c2_eplet_asso.txt")
c2_result <- read_table("c2_eplet_asso_withouthHLAms.txt")
head(c2_result)
min(c2_result$P)



head(c2_result)

grep("DR",c2_result$ID)
grepl("DR",c2_result$ID)
c2_result %>% as.data.frame() %>% mutate('log10P' = log10(as.numeric(P))) %>% #max(log10P)
  mutate('type' = ifelse(grepl("DR",ID),"DR",ifelse(grepl("DQ",ID),"DQ",ifelse(grepl("DP",ID),"DP","Total")))) -> c2_p


c2_result %>% as.data.frame() %>% mutate('log10P' = log10(as.numeric(P))) %>% #max(log10P)
  mutate('type' = ifelse(grepl("DR",ID),"DR",ifelse(grepl("DQ",ID),"DQ",ifelse(grepl("DP",ID),"DP","Total")))) %>%
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
  mutate('log10P' = log10(as.numeric(P))) %>% #max(log10P)
  filter(log10P > -1.5) %>%
  mutate('Gene' = ifelse(grepl("DR",ID),"DR",ifelse(grepl("DQ",ID) | grepl("QA",ID),"DQ",ifelse(grepl("DP",ID)|grepl("PA",ID),"DP","Total")))) %>% #head()
  mutate('type' = ifelse(grepl("Ab",ID),"Ab",ifelse(grepl("Ot",ID) | grepl("ot",ID) ,"Other","E"))) %>% #count(type)
  #filter(type == "E")
  ggplot(aes(x=fct_inorder(ID),y=-log10P,color=Gene,shape=type)) +
  geom_point(size = 2) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom")


c2_result %>% filter(!ID %in% c("All_AbV","All_Eplet","Nr_allDRB","Nr_AbDRB","Nr_otDRB","Nr_allDQB",
                                "Nr_AbDQB","Nr_otDQB","Nr_allDQA","Nr_AbDQA","Nr_otDQA","Nr_allDPB","Nr_AbDPB","Nr_otDPB",
                                "Nr_allDPA","Nr_AbDPA","Nr_otDPA")) %>% #head()
  mutate('log10P' = log10(as.numeric(P))) %>% #head()
  #filter(log10P > -1.5) %>%
  mutate('Gene' = ifelse(grepl("DR",ID),"DR",ifelse(grepl("DQ",ID) | grepl("QA",ID),"DQ",ifelse(grepl("DP",ID)|grepl("PA",ID),"DP","Total")))) %>% #head()
  mutate('type' = ifelse(grepl("Ab",ID),"Ab",ifelse(grepl("Ot",ID) | grepl("ot",ID) ,"Other","E"))) -> c2_result_v2

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

head(c2_result_v2)

head(c2_result)
c2_result_v2 %>% filter(P < 0.06) -> a

c2_result %>% filter(ID %in% c("All_AbV","All_Eplet","Nr_allDRB","Nr_allDQB",
                               "Nr_allDQA","Nr_allDPB","Nr_allDPA")) -> a
  ##filter(grepl(c("All","all"),ID))

c2_result %>% filter(ID %in% c("All_AbV","All_Eplet","Nr_allDRB","Nr_AbDRB","Nr_otDRB","Nr_allDQB",
                                "Nr_AbDQB","Nr_otDQB","Nr_allDQA","Nr_AbDQA","Nr_otDQA","Nr_allDPB","Nr_AbDPB","Nr_otDPB",
                                "Nr_allDPA","Nr_AbDPA","Nr_otDPA")) %>% #head()
  filter(P < 0.1) -> a
  
  
c2 %>% count(AbDR_30C)

c2_result %>% 
  mutate('log10P' = log10(as.numeric(P))) %>% filter(-log10P > 1.5) -> a


#####  EMS eplet weighted sum

c1 <- read.table("../Association/c1_eplet_asso_withouthHLAms.txt",header = T)
c2 <- read.table("../Association/c2_eplet_asso_withouthHLAms.txt",header = T)


head(c1)
head(c2)
head(c1_data)
head(c2_data)



c1_aa <- c1[4:nrow(c1),]
c2_aa <- c2[18:nrow(c2),]
head(c1_aa)
head(c2_aa)

c1_aa %>% mutate('Effect' = Estimate) %>%
  mutate('Effect(t1)' = ifelse(P<0.05,Estimate,0)) %>%
  mutate('Effect(t2)' = ifelse(P<0.01,Estimate,0)) %>% 
  select(ID,Effect,`Effect(t1)`,`Effect(t2)`)-> c1_aa

c2_aa %>% mutate('Effect' = Estimate) %>%
  mutate('Effect(t1)' = ifelse(P<0.05,Estimate,0)) %>%
  mutate('Effect(t2)' = ifelse(P<0.01,Estimate,0)) %>%
  select(ID,Effect,`Effect(t1)`,`Effect(t2)`)-> c2_aa

head(c1_aa)
head(c2_aa)


c1_data %>% select(bCODE, c1_aa$ID) -> c1_aa_data
c2_data %>% select(bCODE, c2_aa$ID) -> c2_aa_data



head(c1_aa)
head(c2_aa)

head(c1_aa)
head(c1_aa_data[1:5])
head(c2_aa_data)


result = NULL
for (i in 1:nrow(c1_aa_data)) {
  result <- rbind(result,c(c1_aa_data$bCODE[i],sum(c1_aa[,2] * c1_aa_data[i,2:ncol(c1_aa_data)]),sum(c1_aa[,3] * c1_aa_data[i,2:ncol(c1_aa_data)]),sum(c1_aa[,4] * c1_aa_data[i,2:ncol(c1_aa_data)])))
}
colnames(result) <- c("bCODE","Effect","Effect_t1","Effect_t2")
c1_result<- result %>% as.data.frame()

result = NULL
for (i in 1:nrow(c1_aa_data)) {
  result <- rbind(result,c(c2_aa_data$bCODE[i],sum(c2_aa[,2] * c2_aa_data[i,2:ncol(c2_aa_data)]),sum(c2_aa[,3] * c2_aa_data[i,2:ncol(c2_aa_data)]),sum(c2_aa[,4] * c2_aa_data[i,2:ncol(c2_aa_data)])))
}
colnames(result) <- c("bCODE","Effect","Effect_t1","Effect_t2")
c2_result<- result %>% as.data.frame()


head(c1_result)
head(c2_result)

#write.table(c1_result,"c1_epletAA_weigthedSum.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(c2_result,"c2_epletAA_weigthedSum.txt",col.names = T,row.names = F,quote = F,sep = "\t")

c1_result <- read.table("c1_epletAA_weigthedSum.txt",header = T)
c2_result <- read.table("c2_epletAA_weigthedSum.txt",header = T)

head(c1_result)
head(c2_result)
colnames(c1_result) <- c("bCODE","Effect","Effect1","Effect2")
colnames(c2_result) <- c("bCODE","Effect","Effect1","Effect2")
head(pheno)
pheno %>% merge(ref,by.x = "KCHIP_ID",by.y="KBA_ID") %>% 
  merge(c1_result) -> c1_data2

pheno %>% merge(ref,by.x = "KCHIP_ID",by.y="KBA_ID") %>% 
  merge(c2_result) -> c2_data2

head(c1_data2)
head(c2_data2)


result <- NULL
for(i in c("Effect","Effect1","Effect2")){
#  print(i)
  #table(c2_data[,i])
  temp <- glm(paste("rej_tot ~ ",i, "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c1_data2, family="binomial")
  result <- rbind(result, c(i,as.vector(summary(temp)$coefficients[2,])))
}

head(result)
colnames(result) <- c("ID","Estimate","SE","Z","P")
c1_result2 <- result %>% as.data.frame()
head(c1_result2)

result <- NULL
for(i in c("Effect","Effect1","Effect2")){
  #table(c2_data[,i])
  temp <- glm(paste("rej_tot ~ ",i, "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_data2, family="binomial")
  result <- rbind(result, c(i,as.vector(summary(temp)$coefficients[2,])))
}
#result <- rbind(result, c(colnames(c2_data)[i],as.vector(summary(temp)$coefficients[2,])))
colnames(result) <- c("ID","Estimate","SE","Z","P")
c2_result2 <- result %>% as.data.frame()

head(c1_result2)
head(c2_result2)
c1_result2$class <- "classI"
c2_result2$class <- "classII"

#write.table(rbind(c1_result2,c2_result2),"c1c2_eplet_weigthSum_association.Result.txt",col.names = T,row.names = F,quote = F,sep = "\t")


head(c1_result)
head(c1_data)
head(c1_result)

result <- NULL
temp <- glm(paste("rej_tot ~ ","Effect", "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_data2, family="binomial")
result <- rbind(result, c(i,as.vector(summary(temp)$coefficients[2,])))

temp <- glm(paste("rej_tot ~ ","Effect1", "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_data2, family="binomial")
result <- rbind(result, c(i,as.vector(summary(temp)$coefficients[2,])))

temp <- glm(paste("rej_tot ~ ","Effect2", "+AGE+SEX+dm+cvd+cmv_igg_reci+hbsag_reci+hcv_ab_reci+ind_atg+D_AGE+D_SEX+dgf",sep=""), data=c2_data2, family="binomial")
result <- rbind(result, c(i,as.vector(summary(temp)$coefficients[2,])))

result

rbind(c1_result2,c2_result2)
#writexl::write_xlsx(c1_data,'test.xlsx')
