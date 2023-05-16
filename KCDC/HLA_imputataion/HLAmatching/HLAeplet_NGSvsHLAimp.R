### HLA eplet about NGS vs
### 1. HLA NGS typing check
### 2. HLA typing vs HLA imputation
library(tidyverse)
library(ggpubr)


##### 20230508 HLA imp result vs NGS result
dr <- read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/QClist/Allogenomics_KR.KD.QCin.pairTable.txt",header = T)
ref <- read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/QClist/Allogenomics_KR.KD.QCin.pairTable.txt",header = T)
head(ref)
ref %>% select(KBA_ID.x,bCODE.x,KBA_ID.y,bCODE.y) -> ref
colnames(ref) <- c("KR","KR_bCODE","KD","KD_bCODE")
head(ref)
head(dr)
dr %>% select(KBA_ID.x,KBA_ID.y) -> dr
colnames(dr) <- c("KR","KD")
head(dr)

ngs <- read.table("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_529sample_IMGT3320convert.txt",header = T)
ngs %>% select(-Sample) -> ngs
head(ngs)
head(dr)


dr %>% merge(ngs,by.x = 'KR',by.y = "ID") %>%
  merge(ngs,by.x = 'KD',by.y = "ID") -> ngs_out

head(ngs_out)
tcol <- grep("NGS_",colnames(ngs_out))
tcol
ngs_out[,tcol] <- lapply(ngs_out[tcol], function(x)
  paste0(str_split_fixed(x,":",3)[,1],":",str_split_fixed(x,":",3)[,2])
)
head(ngs_out)
ngs_out$x <- "x"
ngs_out$x1 <- "x"
ngs_out$x2 <- "x"
ngs_out$x3 <- "x"

ngs_out %>%
  select(KR,NGS_A.1.x,NGS_A.2.x,NGS_B.1.x,NGS_B.2.x,NGS_C.1.x,NGS_C.2.x,KD,NGS_A.1.y,NGS_A.2.y,NGS_B.1.y,NGS_B.2.y,NGS_C.1.y,NGS_C.2.y) -> c1
#writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class1.xlsx",col_names = F,) -> c1
ngs_out %>%
  select(KR,NGS_DRB1.1.x,NGS_DRB1.2.x,x,x1,NGS_DQB1.1.x,NGS_DQB1.2.x,NGS_DQA1.1.x,NGS_DQA1.2.x,NGS_DPB1.1.x,NGS_DPB1.2.x,NGS_DPA1.1.x,NGS_DPA1.2.x,KD,NGS_DRB1.1.y,NGS_DRB1.2.y,x2,x3,NGS_DQB1.1.y,NGS_DQB1.2.y,NGS_DQA1.1.y,NGS_DQA1.2.y,NGS_DPB1.1.y,NGS_DPB1.2.y,NGS_DPA1.1.y,NGS_DPA1.2.y) -> c2
#writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class2.xlsx",col_names = F) -> c2
dataset_names <- list("classI" = c1, "classII" = c2)
writexl::write_xlsx(dataset_names,"test_mtsheet.xlsx")




library(readxlsb)
setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

df <- read_xlsb("NGS_eplet/ABC_Eplet_Matching_4.0_HLAtyping.xlsb",sheet = "Check HLA")
head(df)

colnames(df)
select(grep("Check",colnames(df)))
df %>% select(grep("Check",colnames(df))) %>% #head()
  filter(!row_number() %in% c(1,2)) %>% #head()
  pivot_longer(1:12,names_to = "check",values_to = "type") %>% #head()
  na.omit() %>% select(type) %>% unique() -> c1

df <- read_xlsb("NGS_eplet/DRDQDP_Eplet_Matching_3.1_HLAtyping.xlsb",sheet = "Check HLA")
#colnames(df)
df %>% select(grep("check",colnames(df))) %>% #head()
  filter(!row_number() %in% c(1,2)) %>% #head()
  pivot_longer(1:24,names_to = "check",values_to = "type") %>% #head()
  na.omit() %>% select(type) %>% #head()
  unique() -> c2


head(c1)
head(c2)

c1$class <- "classI"
c2$class <- "classII"

c <- rbind(c1,c2)
head(c)

c %>% select(class,type) %>% arrange(class,type) %>% na.omit() %>% head()
  #write_xlsx("NGS_eplet/Check.HLAtypelist.foreplet.xlsx")

c %>% select(class,type) %>% arrange(class,type) %>% na.omit() -> c_e 



### 2. HLA typing vs HLA imputation
### after HLA type check 20230510
library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")


ngs <- read.table("~/Desktop/KCDC/HLAimputation/HLAtyping_Final_20211126/IMGT3320/HLA.typing.Final.result_529sample_IMGT3320convert.txt",header = T)
ngs %>% select(-Sample) -> ngs
head(ngs)
head(dr)
head(ref)

ref %>% merge(ngs,by.x = 'KR',by.y = "ID") %>%
  merge(ngs,by.x = 'KD',by.y = "ID") -> ngs_out

head(ngs_out)
tcol <- grep("NGS_",colnames(ngs_out))
tcol
ngs_out[,tcol] <- lapply(ngs_out[tcol], function(x)
  paste0(str_split_fixed(x,":",3)[,1],":",str_split_fixed(x,":",3)[,2])
)

check_list <- read_excel("NGS_eplet/Check.HLAtypelist (version 2)_????????????????????????????????????????????????.xlsx")
head(check_list)
check_list %>% filter(type != "DQA1*02:01") %>% select(type,...3) -> check_list

'
ngs_out %>% mutate_all(funs(str_replace(.,"C\\*08:41","C\\*08:01"))) %>% 
  mutate_all(funs(str_replace(.,"DRB1\\*14:141","DRB1\\*14:03"))) %>%
  mutate_all(funs(str_replace(.,"DRB1\\*12:17","DRB1\\*12:01"))) %>%
  mutate_all(funs(str_replace(.,"DQA1\\*05:06","DQA1\\*05:03"))) %>%
  mutate_all(funs(str_replace(.,"DQA1\\*05:08","DQA1\\*05:05"))) -> hlaimp
'
colnames(check_list)[2] <- "to"
head(check_list)

for (i in 1:nrow(check_list)) {
  #ngs_out <- ngs_out %>% mutate_all(funs(str_replace(.,check_list$type[i],check_list$to[i])))
  a <- str_split_fixed(check_list$type[i],"\\*",2)
  b <- str_split_fixed(check_list$to[i],"\\*",2)
  ngs_out <- ngs_out %>% mutate_all(funs(str_replace(.,paste0(a[,1],"\\*",a[,2]),paste0(b[,1],"\\*",b[,2]))))
  #ngs_out[ngs_out == paste0(a[,1],"\\*",a[,2])] <-  paste0(b[,1],"\\*",b[,2])
}

ngs_out <- ngs_out %>% mutate_all(funs(str_replace(.,"x\\*","x"))) %>%
  mutate_all(funs(str_replace(.,"^0:$","x"))) 

head(ngs_out)
ngs_out$x <- "x"
ngs_out$x1 <- "x"
ngs_out$x2 <- "x"
ngs_out$x3 <- "x"

ngs_out %>%
  select(KR_bCODE,NGS_A.1.x,NGS_A.2.x,NGS_B.1.x,NGS_B.2.x,NGS_C.1.x,NGS_C.2.x,KD_bCODE,NGS_A.1.y,NGS_A.2.y,NGS_B.1.y,NGS_B.2.y,NGS_C.1.y,NGS_C.2.y) -> c1
#writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class1.xlsx",col_names = F,) -> c1
ngs_out %>%
  select(KR_bCODE,NGS_DRB1.1.x,NGS_DRB1.2.x,x,x1,NGS_DQB1.1.x,NGS_DQB1.2.x,NGS_DQA1.1.x,NGS_DQA1.2.x,NGS_DPB1.1.x,NGS_DPB1.2.x,NGS_DPA1.1.x,NGS_DPA1.2.x,KD_bCODE,NGS_DRB1.1.y,NGS_DRB1.2.y,x2,x3,NGS_DQB1.1.y,NGS_DQB1.2.y,NGS_DQA1.1.y,NGS_DQA1.2.y,NGS_DPB1.1.y,NGS_DPB1.2.y,NGS_DPA1.1.y,NGS_DPA1.2.y) -> c2
#writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class2.xlsx",col_names = F) -> c2

dataset_names <- list("classI" = c1, "classII" = c2)
writexl::write_xlsx(dataset_names,"NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.afterHLAtypecheck.bCODE.xlsx",col_names = F)

#### after input HLA matchmaker

c1_ngs <- read_xlsb("NGS_eplet/ABC_Eplet_Matching_4.0_HLAtyping.xlsb",sheet = 4,skip = 4)
c1_ngs[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-") %>% tail()
c1_ngs <- c1_ngs[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")


c2_ngs <- read_xlsb("NGS_eplet/DRDQDP_Eplet_Matching_3.1_HLAtyping.xlsb",sheet = 4,skip = 4)
c2_ngs <- c2_ngs[,c(1:54)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")

head(c1_ngs)
head(c2_ngs)

#class I
#KR	Rec1stA	Rec2ndA	Rec1stB	Rec2ndB	Rec1stC	Rec2ndC	KD	D1stA	D2ndA	D1stB	D2ndB	D1stC	D2ndC	Outcome	All Eplets	AbVEp	OthEp	AbVer Eplets	Other Eplets
colnames(c1_ngs) <- c("RecInfo","Rec1stA","Rec2ndA","Rec1stB","Rec2ndB","Rec1stC","Rec2ndC","DonInfo","D1stA","D2ndA","D1stB","D2ndB","D1stC","D2ndC","Outcome","All_Eplets","AbVEp","OthEp","AbVer_Eplets","Other_Eplets")
#class II
colnames(c2_ngs) <- c("RecInfo","Rec1stDRB","Rec2ndDRB","Rec1stDRW","Rec2ndDRW","Rec1stDQB","Rec2ndDQB","Rec1stDQA",
                  "Rec2ndDQA","Rec1stDPB","Rec2ndDPB","Rec1stDPA","Rec2ndDPA","DonInfo","D1stDRB","D2ndDRB","D1stDRW",
                  "D2ndDRW","D1stDQB","D2ndDQB","D1stDQA","D2ndDQA","D1stDPB","D2ndDPB","D1stDPA","D2ndDPA",	
                  "About_Outcome","All_ABV_ClassII","Total_eps","allDRB","AbDRB","otDRB","allDQB","AbDQB","otDQB",
                  "allDQA","AbDQA","otDQA","allDPB","AbDPB","otDPB","allDPA","AbDPA","otDPA","AbDRB","OtDRB","AbDQB",	
                  "OtDQB","AbDQA","OtDQA","AbDPB","OtDPB","AbDPA","OtDPA")

head(c2_ngs)


##
c1_imp <- read_xlsx("Association/c1_forWAS.xlsx")
c2_imp <- read_xlsx("Association/c2_forWAS.xlsx")
head(c1_imp)
colnames(c1_imp)
colnames(c2_imp)


c1_ngs %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>% mutate("type" = "NGS") %>% head()
c1_imp %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>% mutate("type" = "HLAimp") %>% filter(RecInfo %in% c1_ngs$RecInfo)

c1_ngs %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>% mutate("type" = "NGS") %>% 
  rbind(c1_imp %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>% mutate("type" = "HLAimp") %>% filter(RecInfo %in% c1_ngs$RecInfo)) %>%
  pivot_longer(3:5,names_to = "theme",values_to = "count") %>% #head()
  pivot_wider(names_from = type,values_from = count) %>%
  ggplot(aes(x=NGS,y=HLAimp,color=theme)) + 
  geom_point() +
  facet_grid(~theme) + 
  theme(legend.position = "none")


c1_ngs %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>% mutate("type" = "NGS") %>% 
  rbind(c1_imp %>% select(RecInfo,DonInfo,All_Eplets,AbVEp,OthEp) %>% mutate("type" = "HLAimp") %>% filter(RecInfo %in% c1_ngs$RecInfo)) %>%
  pivot_longer(3:5,names_to = "theme",values_to = "count") %>% #head()
  pivot_wider(names_from = type,values_from = count) %>%
  mutate(theme = factor(theme,levels = c("All_Eplets","AbVEp","OthEp"))) %>%
  #ggplot(aes(x=NGS,y=HLAimp,color=theme)) + 
  ggscatter(.,x='NGS',y='HLAimp',color='theme',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            facet.by = "theme",
            xlab = "NGS based HLA typing",
            ylab = "HLA imputation",
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) + 
  theme(legend.position = "none",
        strip.text.x = element_text(size = 13,face = "bold"))

head(c2_ngs)
head(c2_imp)
colnames(c2_imp)[28:44]
colnames(c2_ngs)[28:44]
c2_imp <- c2_imp %>% select(RecInfo,DonInfo,28:44) %>% mutate("type" = "HLAimp") %>% filter(RecInfo %in% c2_ngs$RecInfo)
colnames(c2_imp) <- str_replace_all(colnames(c2_imp),"Nr_","")
colnames(c2_imp)

c2_ngs %>% select(RecInfo,DonInfo,28:44) %>% mutate("type" = "NGS") %>% #head()
  rbind(c2_imp) %>%
  pivot_longer(3:19,names_to = "theme",values_to = "count") %>% #head()
  pivot_wider(names_from = type,values_from = count) %>% 
  filter(theme %in% c("Total_eps","All_ABV_ClassII","allDPB","allDPA","allDQB","allDQA","allDRB")) %>% #count(theme)
  mutate(theme = factor(theme,levels = c("Total_eps","All_ABV_ClassII","allDPB","allDPA","allDQB","allDQA","allDRB"))) %>%
  #ggplot(aes(x=NGS,y=HLAimp,color=theme)) + 
  filter(theme %in% c("Total_eps","All_ABV_ClassII")) %>%
  ggscatter(.,x='NGS',y='HLAimp',color='theme',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            facet.by = "theme",
            xlab = "NGS based HLA typing",
            ylab = "HLA imputation",nrow =1,
            scales='free',
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "none",
        strip.text.x = element_text(size = 13,face = "bold"))


c2_ngs %>% select(RecInfo,DonInfo,28:44) %>% mutate("type" = "NGS") %>% #head()
  rbind(c2_imp) %>%
  pivot_longer(3:19,names_to = "theme",values_to = "count") %>% #head()
  pivot_wider(names_from = type,values_from = count) %>% 
  filter(!theme %in% c("Total_eps","All_ABV_ClassII","allDRB","allDQB","allDQA","allDPB","allDPA")) %>% #count(theme)
  #mutate(theme = factor(theme,levels = c("Total_eps","All_ABV_ClassII","allDRB","allDQB","allDQA","allDPB","allDPA"))) %>%
  #ggplot(aes(x=NGS,y=HLAimp,color=theme)) + 
  ggscatter(.,x='NGS',y='HLAimp',color='theme',
            add = "reg.line",
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            facet.by = "theme",
            xlab = "NGS based HLA typing",
            ylab = "HLA imputation",nrow =2,scales='free',
            #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
            cor.coeff.args = list(method = "pearson", label.sep = "\n")) + 
  theme(legend.position = "none",
        strip.text.x = element_text(size = 12,face = "bold"))
