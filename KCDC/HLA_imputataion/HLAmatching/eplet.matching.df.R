## eplet matching data frame 
library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

hlaimp <- read.table("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/test/cookHLA_han/KRKD.HLAimp_cookHLA.HANref.hg18.MHC.HLA_IMPUTATION_OUT.Nomencleaner_withheader.chped",header = T)
head(hlaimp)

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
head(hlaimp)

#class1
'
RecInfo	Rec1stA	Rec2ndA	Rec1stB	Rec2ndB	Rec1stC	Rec2ndC	DonorInfo	Don1stA	Don2ndA	Don1stB	Don2ndB	Don1stC	Don2ndC
'
'
ID	1stDRB	2ndDRB	1stDRW	2ndDRW	1stDQB	2ndDQB	1stDQA	2ndDQA	1stDPB	2ndDPB	1stDPA	2ndDPA	
ID	1stDRB	2ndDRB	1stDRW	2ndDRW	1stDQB	2ndDQB	1stDQA	2ndDQA	1stDPB	2ndDPB	1stDPA	2ndDPA
'


### replace for eplet metching
#df %>% pivot_longer(cols = 3:ncol(df),names_to = "Gene",values_to = 'type') %>%
  filter(type %in% c("C*08:41"," DRB1*14:141","DRB1*12:17","DQA1*05:06","DQA1*05:08",
                     "C*08:01","DRB1*14:03","DRB1*12:01","DQA1*05:03","DQA1*05:05"))
  
  
hlaimp %>% replace

hlaimp  %>% mutate_all(funs(str_repl))
hlaimp %>% filter(HLA_C.2 == 'C*08:41')


hlaimp %>% mutate_all(funs(str_replace(.,"C\\*08:41","C\\*08:01"))) %>% #filter(IID %in% c('NIH19KT0031','NIH19KT6205',"NIH19KT6315","NIH19KT6426"))
  mutate_all(funs(str_replace(.,"DRB1\\*14:141","DRB1\\*14:03"))) %>%
  mutate_all(funs(str_replace(.,"DRB1\\*12:17","DRB1\\*12:01"))) %>%
  mutate_all(funs(str_replace(.,"DQA1\\*05:06","DQA1\\*05:03"))) %>%
  mutate_all(funs(str_replace(.,"DQA1\\*05:08","DQA1\\*05:05"))) -> hlaimp

###

hlaimp %>% select(IID,HLA_A.1,HLA_A.2,HLA_B.1,HLA_B.2,HLA_C.1,HLA_C.2) -> hla_c1
hlaimp %>% select(IID,HLA_DRB1.1,HLA_DRB1.2,HLA_DQB1.1,HLA_DQB1.2,HLA_DQA1.1,HLA_DQA1.2,HLA_DPB1.1,HLA_DPB1.2,HLA_DPA1.1,HLA_DPA1.2) -> hla_c2
head(hla_c1)
head(hla_c2)
#colnames(hla_c1)[1] <- "KR"
#colnames(hla_c2)[1] <- "KR"

head(dr)

#dr %>% left_join(hla_c1) -> hla_c1_epi

##### 1000개 안나누고 저장
dr %>% merge(hla_c1,by.x = "KR",by.y = "IID") %>% merge(hla_c1,by.x = "KD",by.y = "IID") %>% 
  select(KR,HLA_A.1.x,HLA_A.2.x,HLA_B.1.x,HLA_B.2.x,HLA_C.1.x,HLA_C.2.x,KD,HLA_A.1.y,HLA_A.2.y,HLA_B.1.y,HLA_B.2.y,HLA_C.1.y,HLA_C.2.y) %>%
  writexl::write_xlsx("KR.KD.HLAimp.foreplet.class1.v2.xlsx",col_names = F)


dr %>% merge(hla_c2,by.x = "KR",by.y = "IID") %>% merge(hla_c2,by.x = "KD",by.y = "IID") %>%  #head()
  mutate("x" = "x","x1" = "x","x2" = "x","x3" = "x") %>% #head()
  select(KR,HLA_DRB1.1.x,HLA_DRB1.2.x,x,x1,HLA_DQB1.1.x,HLA_DQB1.2.x,HLA_DQA1.1.x,HLA_DQA1.2.x,HLA_DPB1.1.x,HLA_DPB1.2.x,HLA_DPA1.1.x,HLA_DPA1.2.x,KD,HLA_DRB1.1.y,HLA_DRB1.2.y,x2,x3,HLA_DQB1.1.y,HLA_DQB1.2.y,HLA_DQA1.1.y,HLA_DQA1.2.y,HLA_DPB1.1.y,HLA_DPB1.2.y,HLA_DPA1.1.y,HLA_DPA1.2.y) %>%
  writexl::write_xlsx("KR.KD.HLAimp.foreplet.class2.v2.xlsx",col_names = F)

##### 1000개 나누고 한 파일에 저장 ## save by bCODE
ref %>% merge(hla_c1,by.x = "KR",by.y = "IID") %>% #dim()
  merge(hla_c1,by.x = "KD",by.y = "IID") %>% #head()
  #select(KR,HLA_A.1.x,HLA_A.2.x,HLA_B.1.x,HLA_B.2.x,HLA_C.1.x,HLA_C.2.x,KD,HLA_A.1.y,HLA_A.2.y,HLA_B.1.y,HLA_B.2.y,HLA_C.1.y,HLA_C.2.y) -> c1
  select(KR_bCODE,HLA_A.1.x,HLA_A.2.x,HLA_B.1.x,HLA_B.2.x,HLA_C.1.x,HLA_C.2.x,KD_bCODE,HLA_A.1.y,HLA_A.2.y,HLA_B.1.y,HLA_B.2.y,HLA_C.1.y,HLA_C.2.y) -> c1

ref %>% merge(hla_c2,by.x = "KR",by.y = "IID") %>% merge(hla_c2,by.x = "KD",by.y = "IID") %>%  #head()
  mutate("x" = "x","x1" = "x","x2" = "x","x3" = "x") %>% #head()
  #select(KR,HLA_DRB1.1.x,HLA_DRB1.2.x,x,x1,HLA_DQB1.1.x,HLA_DQB1.2.x,HLA_DQA1.1.x,HLA_DQA1.2.x,HLA_DPB1.1.x,HLA_DPB1.2.x,HLA_DPA1.1.x,HLA_DPA1.2.x,KD,HLA_DRB1.1.y,HLA_DRB1.2.y,x2,x3,HLA_DQB1.1.y,HLA_DQB1.2.y,HLA_DQA1.1.y,HLA_DQA1.2.y,HLA_DPB1.1.y,HLA_DPB1.2.y,HLA_DPA1.1.y,HLA_DPA1.2.y) -> c2
  select(KR_bCODE,HLA_DRB1.1.x,HLA_DRB1.2.x,x,x1,HLA_DQB1.1.x,HLA_DQB1.2.x,HLA_DQA1.1.x,HLA_DQA1.2.x,HLA_DPB1.1.x,HLA_DPB1.2.x,HLA_DPA1.1.x,HLA_DPA1.2.x,KD_bCODE,HLA_DRB1.1.y,HLA_DRB1.2.y,x2,x3,HLA_DQB1.1.y,HLA_DQB1.2.y,HLA_DQA1.1.y,HLA_DQA1.2.y,HLA_DPB1.1.y,HLA_DPB1.2.y,HLA_DPA1.1.y,HLA_DPA1.2.y) -> c2


dataset_names <- list("classI_g1" = c1[1:1000,], "classI_g2" = c1[1001:nrow(c1),],"classII_g1" = c2[1:1000,], "classII_g2" = c2[1001:nrow(c2),])
writexl::write_xlsx(dataset_names,"KR.KD.HLAimp.foreplet.v2_bCODE.xlsx",col_names = F)


c1[1001:nrow(c1),] %>% dim()






######
library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230420_eplet/")

c1_1 <- read_xlsb("ABC_Eplet_Matching_4.0_g1.xlsb",sheet = 4,skip = 4)
c1_2 <- read_xlsb("ABC_Eplet_Matching_4.0_g2.xlsb",sheet = 4,skip = 4)
colnames(c1_1)
#df %>% select(grep("column",colnames(df)),X,X.1,X.2,X.3,X.4,X.5,X.6) 
c1_1[,c(1:20)] %>% #tail()
  filter(X != "x") %>% filter(!is.na(X)) %>% #tail()  
  filter(X != "") %>%  head()
  dim()

c1_2[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-") %>% tail()
  
c1_1 <- c1_1[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")
c1_2 <- c1_2[,c(1:20)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")
tail(c1_1)
head(c1_1)
c1 <- rbind(c1_1,c1_2)

c2_1 <- read_xlsb("DRDQDP_Eplet_Matching_3.1_g1.xlsb",sheet = 4,skip = 4)
c2_2 <- read_xlsb("DRDQDP_Eplet_Matching_3.1_g2.xlsb",sheet = 4,skip = 4)
colnames(c2_2)

head(c2_1)
head(c2_2)

  
c2_1 <- c2_1[,c(1:55)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")
c2_2 <- c2_2[,c(1:55)] %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")

head(c2_1)
head(c2_2)


c2 <- rbind(c2_1,c2_2)
head(c2)
#class I
#KR	Rec1stA	Rec2ndA	Rec1stB	Rec2ndB	Rec1stC	Rec2ndC	KD	D1stA	D2ndA	D1stB	D2ndB	D1stC	D2ndC	Outcome	All Eplets	AbVEp	OthEp	AbVer Eplets	Other Eplets
colnames(c1) <- c("KR","Rec1stA","Rec2ndA","Rec1stB","Rec2ndB","Rec1stC","Rec2ndC","KD","D1stA","D2ndA","D1stB","D2ndB","D1stC","D2ndC","Outcome","All Eplets","AbVEp","OthEp","AbVer Eplets","Other Eplets")
#class II
colnames(c2) <- c("KR","Rec1stDRB","Rec2ndDRB","Rec1stDRW","Rec2ndDRW","Rec1stDQB","Rec2ndDQB","Rec1stDQA",
                  "Rec2ndDQA","Rec1stDPB","Rec2ndDPB","Rec1stDPA","Rec2ndDPA","KD","D1stDRB","D2ndDRB","D1stDRW",
                  "D2ndDRW","D1stDQB","D2ndDQB","D1stDQA","D2ndDQA","D1stDPB","D2ndDPB","D1stDPA","D2ndDPA",	
                  "Outcome","ALL AbV","ClassII","total","allDRB","AbDRB","otDRB","allDQB","AbDQB","otDQB",
                  "allDQA","AbDQA","otDQA","allDPB","AbDPB","otDPB","allDPA","AbDPA","otDPA","AbDRB","OtDRB","AbDQB",	
                  "OtDQB","AbDQA","OtDQA","AbDPB","OtDPB","AbDPA","OtDPA")
#KR	Rec1stDRB	Rec2ndDRB	Rec1stDRW	Rec2ndDRW	Rec1stDQB	Rec2ndDQB	Rec1stDQA	Rec2ndDQA	Rec1stDPB	Rec2ndDPB	Rec1stDPA	Rec2ndDPA	KD	D1stDRB	D2ndDRB	D1stDRW	D2ndDRW	D1stDQB	D2ndDQB	D1stDQA	D2ndDQA	D1stDPB	D2ndDPB	D1stDPA	D2ndDPA	Outcome	ALL AbV ClassII	total	allDRB	AbDRB	otDRB	allDQB	AbDQB	otDQB	allDQA	AbDQA	otDQA	allDPB	AbDPB	otDPB	allDPA	AbDPA	otDPA	AbDRB	OtDRB	AbDQB	OtDQB	AbDQA	OtDQA	AbDPB	OtDPB	AbDPA	OtDPA

#writexl::write_xlsx(c1,"KR.KD.HLAimp.foreplet.class1_info.xlsx",col_names = F)
#writexl::write_xlsx(c2,"KR.KD.HLAimp.foreplet.class2_info.xlsx",col_names = F)


c1 <- read_excel("20230420_eplet/KR.KD.HLAimp.foreplet.classI_ABC_info.xlsx")
c2 <- read_excel("20230420_eplet/KR.KD.HLAimp.foreplet.classII_DRDQDP_info.xlsx")

head(c2)
head(c1)
colnames(c1)
c1 %>% select(KR,KD,"All Eplets",AbVEp,OthEp) %>%
  pivot_longer(cols = 3:5,names_to = "type",values_to = "score") %>%
  ggplot(aes(x=score,fill=type))+
  geom_histogram() +
  facet_grid(~type) +
  theme(legend.position = 'none')

colnames(c2)
c2 %>% select(KR,KD,28:44) %>% #dim()
  pivot_longer(cols = 3:4,names_to = "type",values_to = "score") %>%
  ggplot(aes(x=score,fill=type))+
  geom_histogram() +
  facet_grid(~type) +
  theme(legend.position = 'none')


c2 %>% select(KR,KD,28:44) %>%# head()#dim()
  pivot_longer(cols = c(5,8,11,14,17),names_to = "type",values_to = "score") %>%
  ggplot(aes(x=score,fill=type))+
  geom_histogram() +
  facet_grid(~type) +
  theme(legend.position = 'none')

### bcode ID
library(tidyverse)
library(stringr)
library(readxl)

setwd("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/")

ref <- read.table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt",header = T)
head(ref)

c1 <- read_excel("20230420_eplet/KR.KD.HLAimp.foreplet.classI_ABC_info.xlsx")
c2 <- read_excel("20230420_eplet/KR.KD.HLAimp.foreplet.classII_DRDQDP_info.xlsx")
head(c1)

c1 %>% select(KR,KD) %>% merge(ref,by.x = "KR",by.y = "KBA_ID") %>% #head()
  merge(ref,by.x = "KD",by.y = "KBA_ID") -> ref1 #%>% #head()
head(c1 %>% select(KR,KD) )
head(ref1)

c2 %>% select(KR,KD) %>% merge(ref,by.x = "KR",by.y = "KBA_ID") %>% #head()
  merge(ref,by.x = "KD",by.y = "KBA_ID") -> ref2 #%>% #head()


#cbind(c1 %>% select(KR,KD),ref1) %>% #head()
  
c1 %>% mutate(KR = ref1$bCODE.x) %>% mutate(KD = ref1$bCODE.y) %>%
  writexl::write_xlsx("20230420_eplet/forKOTRY/KR.KD.HLAimp.foreplet.classI_ABC_info_20230421.xlsx")

c2 %>% mutate(KR = ref2$bCODE.x) %>% mutate(KD = ref2$bCODE.y) %>%
  writexl::write_xlsx("20230420_eplet/forKOTRY/KR.KD.HLAimp.foreplet.classII_DRDQDP_info_20230421.xlsx")


### original eplet sheet

head(c1_1)
head(c1_2)
head(c2_1)
head(c2_2)

c1_1 %>% select(X,X.1) %>% merge(ref,by.x = "X",by.y = "KBA_ID") %>% #head()
  merge(ref,by.x = "X.1",by.y = "KBA_ID") %>% #head()
  select(X,bCODE.x,X.1,bCODE.y) -> m_c1_1
  #writexl::write_xlsx("20230420_eplet/forKOTRY/NIH.bcode.MATCHING.sheet.xlsx",)
c1_2 %>% select(X,X.1) %>% merge(ref,by.x = "X",by.y = "KBA_ID") %>% #head()
  merge(ref,by.x = "X.1",by.y = "KBA_ID") %>% #head()
  select(X,bCODE.x,X.1,bCODE.y) -> m_c1_2

c2_1 %>% select(X,X.13) %>% merge(ref,by.x = "X",by.y = "KBA_ID") %>% #head()
  merge(ref,by.x = "X.13",by.y = "KBA_ID") %>% #head()
  select(X,bCODE.x,X.13,bCODE.y) -> m_c2_1
#writexl::write_xlsx("20230420_eplet/forKOTRY/NIH.bcode.MATCHING.sheet.xlsx",)
c2_2 %>% select(X,X.13) %>% merge(ref,by.x = "X",by.y = "KBA_ID") %>% #head()
  merge(ref,by.x = "X.13",by.y = "KBA_ID") %>% #head()
  select(X,bCODE.x,X.13,bCODE.y) -> m_c2_2

dataset_names <- list('c1_g1' = m_c1_1,'c1_g2' = m_c1_2,'c2_g1' = m_c2_1,'c2_g2' = m_c2_2)

#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(dataset_names, file = 'forKOTRY/NIH.bcode.MATHCING.sheet.xlsx') 



#c1 %>% select(KR,KD,`AbVer Eplets`,`Other Eplets`) -> test
#head(test)  


##### 20230508 HLA imp result vs NGS result

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
  
#KR,NGS_A.1.x,NGS_A.2.x,NGS_B.1.x,NGS_B.2.x,NGS_C.1.x,NGS_C.2.x,KD,NGS_A.1.y,NGS_A.2.y,NGS_B.1.y,NGS_B.2.y,NGS_C.1.y,NGS_C.2.y
#KR,NGS_DRB1.1.x,NGS_DRB1.2.x,x,x1,NGS_DQB1.1.x,NGS_DQB1.2.x,NGS_DQA1.1.x,NGS_DQA1.2.x,NGS_DPB1.1.x,NGS_DPB1.2.x,NGS_DPA1.1.x,NGS_DPA1.2.x,KD,NGS_DRB1.1.y,NGS_DRB1.2.y,x2,x3,NGS_DQB1.1.y,NGS_DQB1.2.y,NGS_DQA1.1.y,NGS_DQA1.2.y,NGS_DPB1.1.y,NGS_DPB1.2.y,NGS_DPA1.1.y,NGS_DPA1.2.y


ngs_out %>%
  select(KR,NGS_A.1.x,NGS_A.2.x,NGS_B.1.x,NGS_B.2.x,NGS_C.1.x,NGS_C.2.x,KD,NGS_A.1.y,NGS_A.2.y,NGS_B.1.y,NGS_B.2.y,NGS_C.1.y,NGS_C.2.y) -> c1
  #writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class1.xlsx",col_names = F,) -> c1



ngs_out %>%
  select(KR,NGS_DRB1.1.x,NGS_DRB1.2.x,x,x1,NGS_DQB1.1.x,NGS_DQB1.2.x,NGS_DQA1.1.x,NGS_DQA1.2.x,NGS_DPB1.1.x,NGS_DPB1.2.x,NGS_DPA1.1.x,NGS_DPA1.2.x,KD,NGS_DRB1.1.y,NGS_DRB1.2.y,x2,x3,NGS_DQB1.1.y,NGS_DQB1.2.y,NGS_DQA1.1.y,NGS_DQA1.2.y,NGS_DPB1.1.y,NGS_DPB1.2.y,NGS_DPA1.1.y,NGS_DPA1.2.y) -> c2
  #writexl::write_xlsx("./NGS_eplet/KR.KD.HLAtyping_ngs.foreplet.class2.xlsx",col_names = F) -> c2


dataset_names <- list("classI" = c1, "classII" = c2)
#dataset_names <- list("classI_g1" = c1[1:1000,],"classI_g2" = c1[100], "classII_g1" = c2,"classII_g2" = c2)

writexl::write_xlsx(dataset_names,"test_mtsheet.xlsx")

