## eplet matching data frame  HLA matchmaker AR new 20251016
library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)
setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/hla_matchmaker")

check_matchmaker <- readxl::read_xlsx("check_HLAtype_check_forHLAmatchmaker_byParkprof.xlsx")

head(check_matchmaker)
check_matchmaker %>% filter(hlamatchmaker_type == 0) %>% mutate(new_type = ifelse(str_detect(value,"N"),"x",...3)) %>%
  select(value,new_type) -> change_type

change_type

hlaimp <- read_table("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/HLA_KHUref_AR/output_fd.txt") %>% mutate(ID = str_split_fixed(ID,"=",2)[,1])
head(hlaimp)

'''
ref <- read_table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt")
final_id <- read.table("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/SIRPa/new_pair/KOTRY.AR.SIRPA.txt",header = T)
head(final_id)
head(ref)
dim(ref)
ref %>% filter(KBA_ID %in% final_id$ID) %>%
  left_join(hlaimp %>% rename(KBA_ID = ID)) %>% select(-KBA_ID) %>%
  writexl::write_xlsx("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/KOTRY.2882sample.HLAimp.preprocessing_20251209.xlsx")
'''

ref <- readxl::read_xlsx("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/KBA_QC/KBA.QC.list.AR_1441pairs.20250325.xlsx")
head(ref)

hlaimp %>% pivot_longer(cols = 2:ncol(hlaimp),names_to = "Gene",values_to = 'type') %>% #head()
  filter(type %in% c("A*02:474","A*02:07","A*24:294Q","A*24:02","C*07:151","C*03:02"))

hlaimp %>% pivot_longer(cols = 2:ncol(hlaimp),names_to = "Gene",values_to = 'type') %>% #head()
  filter(type %in% c("A*02:474"))
                  

head(check_matchmaker)
check_matchmaker %>% filter(value %in% c("C*08:41"," DRB1*14:141","DRB1*12:17","DQA1*05:06","DQA1*05:08",
                   "C*08:01","DRB1*14:03","DRB1*12:01","DQA1*05:03","DQA1*05:05"))



for (i in seq_len(nrow(change_type))) {
  hlaimp <- hlaimp %>%
    mutate(across(everything(), ~ str_replace_all(., fixed(change_type$value[i]), change_type$new_type[i])))
}
hlaimp
hlaimp %>% pivot_longer(cols = 2:ncol(hlaimp),names_to = "Gene",values_to = 'type') %>% #head()
  filter(type %in% c("A*02:474"))

hlaimp %>% select(ID,HLA_A.1,HLA_A.2,HLA_B.1,HLA_B.2,HLA_C.1,HLA_C.2) -> hla_c1
hlaimp %>% select(ID,HLA_DRB1.1,HLA_DRB1.2,HLA_DQB1.1,HLA_DQB1.2,HLA_DQA1.1,HLA_DQA1.2,HLA_DPB1.1,HLA_DPB1.2,HLA_DPA1.1,HLA_DPA1.2) -> hla_c2
head(hla_c1)
head(hla_c2)



##### 1000�� ������ �� ���Ͽ� ���� ## save by bCODE


head(ref)

ref %>% merge(hla_c1,by.x = "KBA_ID.KR",by.y = "ID") %>% #dim()
  merge(hla_c1,by.x = "KBA_ID.KD",by.y = "ID") %>% #head()
  #select(KR,HLA_A.1.x,HLA_A.2.x,HLA_B.1.x,HLA_B.2.x,HLA_C.1.x,HLA_C.2.x,KD,HLA_A.1.y,HLA_A.2.y,HLA_B.1.y,HLA_B.2.y,HLA_C.1.y,HLA_C.2.y) -> c1
  select(bCODE.KR,HLA_A.1.x,HLA_A.2.x,HLA_B.1.x,HLA_B.2.x,HLA_C.1.x,HLA_C.2.x,bCODE.KD,HLA_A.1.y,HLA_A.2.y,HLA_B.1.y,HLA_B.2.y,HLA_C.1.y,HLA_C.2.y) -> c1

ref %>% merge(hla_c2,by.x = "KBA_ID.KR",by.y = "ID") %>% merge(hla_c2,by.x = "KBA_ID.KD",by.y = "ID") %>%  #head()
  mutate("x" = "x","x1" = "x","x2" = "x","x3" = "x") %>% #head()
  #select(KR,HLA_DRB1.1.x,HLA_DRB1.2.x,x,x1,HLA_DQB1.1.x,HLA_DQB1.2.x,HLA_DQA1.1.x,HLA_DQA1.2.x,HLA_DPB1.1.x,HLA_DPB1.2.x,HLA_DPA1.1.x,HLA_DPA1.2.x,KD,HLA_DRB1.1.y,HLA_DRB1.2.y,x2,x3,HLA_DQB1.1.y,HLA_DQB1.2.y,HLA_DQA1.1.y,HLA_DQA1.2.y,HLA_DPB1.1.y,HLA_DPB1.2.y,HLA_DPA1.1.y,HLA_DPA1.2.y) -> c2
  select(bCODE.KR,HLA_DRB1.1.x,HLA_DRB1.2.x,x,x1,HLA_DQB1.1.x,HLA_DQB1.2.x,HLA_DQA1.1.x,HLA_DQA1.2.x,HLA_DPB1.1.x,HLA_DPB1.2.x,HLA_DPA1.1.x,HLA_DPA1.2.x,bCODE.KD,HLA_DRB1.1.y,HLA_DRB1.2.y,x2,x3,HLA_DQB1.1.y,HLA_DQB1.2.y,HLA_DQA1.1.y,HLA_DQA1.2.y,HLA_DPB1.1.y,HLA_DPB1.2.y,HLA_DPA1.1.y,HLA_DPA1.2.y) -> c2


dataset_names <- list("classI_g1" = c1[1:1000,], "classI_g2" = c1[1001:nrow(c1),],"classII_g1" = c2[1:1000,], "classII_g2" = c2[1001:nrow(c2),])
writexl::write_xlsx(dataset_names,"KR.KD.HLAimp.foreplet.v2_bCODE.xlsx",col_names = F)




######
library(tidyverse)
library(stringr)
library(readxl)
library(readxlsb)

setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/hla_matchmaker")

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

#writexl::write_xlsx(c1,"KR.KD.HLAimp.foreplet.class1_info.xlsx",col_names = T)
#writexl::write_xlsx(c2,"KR.KD.HLAimp.foreplet.class2_info.xlsx",col_names = T)
#head(c1)

c1 <- read_excel("KR.KD.HLAimp.foreplet.class1_info.xlsx")
c2 <- read_excel("KR.KD.HLAimp.foreplet.class2_info.xlsx")
#head(c2)
#head(c1)




#c1 <- read_excel("KR.KD.HLAimp.foreplet.classI_ABC_info.xlsx")
#c2 <- read_excel("KR.KD.HLAimp.foreplet.classII_DRDQDP_info.xlsx")

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
  facet_grid(~type)
  
  
########## with eplet

setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/hla_matchmaker")

c1_1 <- read_xlsb("ABC_Eplet_Matching_4.0_g1.xlsb",sheet = 4,skip = 4)
c1_2 <- read_xlsb("ABC_Eplet_Matching_4.0_g2.xlsb",sheet = 4,skip = 4)
colnames(c1_1)
#df %>% select(grep("column",colnames(df)),X,X.1,X.2,X.3,X.4,X.5,X.6) 
c1_1[,c(1:20)] %>% #tail()
  filter(X != "x") %>% filter(!is.na(X)) %>% #tail()  
  filter(X != "") %>%  head()


head(c1_1)
head(c1_2)

c1_1 <- c1_1 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "")
c1_2 <- c1_2 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-")


c1_hh <- read_excel("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230508_epletv2/c1_header.xlsx",sheet = 2)
c1_hh
colnames(c1_hh)
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
head(c1_forWAS)


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
# �̰� �Ŀ� ���۾�
'

c2_hh <- read_excel("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230508_epletv2/c2_header.xlsx",sheet = 4)
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


writexl::write_xlsx(c1_forWAS,"c1_forWAS.xlsx")
writexl::write_xlsx(c2_forWAS,"c2_forWAS.xlsx")


###
head(c2_forWAS)
head(c1_forWAS)

####  Single molecule Eplet MM 

c2_1 <- read_xlsx("KR.KD.HLAimp.foreplet.v2_bCODE.xlsx",sheet = 3,col_names = F)
c2_2 <- read_xlsx("KR.KD.HLAimp.foreplet.v2_bCODE.xlsx",sheet = 4,col_names = F)

head(c2_1)
head(c2_2)

c2 <- rbind(c2_1,c2_2)
head(c2)
colnames(c2) <- c("KR","Rec1stDRB","Rec2ndDRB","Rec1stDRW","Rec2ndDRW","Rec1stDQB","Rec2ndDQB","Rec1stDQA",
                    "Rec2ndDQA","Rec1stDPB","Rec2ndDPB","Rec1stDPA","Rec2ndDPA","KD","D1stDRB","D2ndDRB","D1stDRW",
                    "D2ndDRW","D1stDQB","D2ndDQB","D1stDQA","D2ndDQA","D1stDPB","D2ndDPB","D1stDPA","D2ndDPA")

#KR_bCODE,HLA_DRB1.1.x,HLA_DRB1.2.x,x,x1,HLA_DQB1.1.x,HLA_DQB1.2.x,HLA_DQA1.1.x,HLA_DQA1.2.x,HLA_DPB1.1.x,HLA_DPB1.2.x,HLA_DPA1.1.x,HLA_DPA1.2.x,KD_bCODE,HLA_DRB1.1.y,HLA_DRB1.2.y,x2,x3,HLA_DQB1.1.y,HLA_DQB1.2.y,HLA_DQA1.1.y,HLA_DQA1.2.y,HLA_DPB1.1.y,HLA_DPB1.2.y,HLA_DPA1.1.y,HLA_DPA1.2.y
colnames(c2)

c2_h1 <- c2[,c(1,2,3,4,4,6,7,4,4,4,4,4,4,14,15,15,4,4,19,19)]
c2_h2 <- c2[,c(1,2,3,4,4,6,7,4,4,4,4,4,4,14,16,16,4,4,20,20)]

head(c2_h1)

dataset_names <- list("classII_g1_h1" = c2_h1[1:1000,], "classII_g2_h1" = c2_h1[1001:nrow(c2_h1),],"classII_g1_h2" = c2_h2[1:1000,], "classII_g2_h2" = c2_h2[1001:nrow(c2_h2),])
writexl::write_xlsx(dataset_names,"SMMM/KR.KD.HLAimp.foreplet.v2_bCODE_singleMM.xlsx",col_names = F)


##### eplet score comparison 
#### after input HLA match maker
c2_g1_h1 <- read_xlsb("SMMM/DRDQDP_Eplet_Matching_3.1_g1_h1_smmm.xlsb",sheet = 4,skip = 4)
c2_g1_h2 <- read_xlsb("SMMM/DRDQDP_Eplet_Matching_3.1_g1_h2_smmm.xlsb",sheet = 4,skip = 4)
c2_g2_h1 <- read_xlsb("SMMM/DRDQDP_Eplet_Matching_3.1_g2_h1_smmm.xlsb",sheet = 4,skip = 4)
c2_g2_h2 <- read_xlsb("SMMM/DRDQDP_Eplet_Matching_3.1_g2_h2_smmm.xlsb",sheet = 4,skip = 4)


head(c2_g1_h1)
head(c2_g1_h2)
head(c2_g2_h1)
head(c2_g2_h2)

c2_g1_h1 <- c2_g1_h1 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "0")
c2_g1_h2 <- c2_g1_h2 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "0")

c2_g2_h1 <- c2_g2_h1 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-") %>% filter(X != "0")
c2_g2_h2 <- c2_g2_h2 %>% filter(X != "x") %>% filter(!is.na(X)) %>% filter(X != "") %>% filter(X != "-") %>% filter(X != "0")

dim(c2_g2_h1)
dim(c2_g2_h2)
dim(c2_g1_h2)
dim(c2_g1_h1)

tail(c2_g1_h2)
tail(c2_g2_h2)
dim(c2_g1_h1)


c2_h1 <- rbind(c2_g1_h1,c2_g2_h1)
c2_h2 <- rbind(c2_g1_h2,c2_g2_h2)

head(c2_h2)
head(c2_h1)

dim(c2_h2)
dim(c2_h1)

'
c2_header <- read_xlsb("DRDQDP_Eplet_Matching_3.1.xlsb",sheet = 4,skip = 1)
#c2_header[1,] %>% t()

c2_header[1,] %>% write_xlsx("c2_header_extra.xlsx")
# �̰� �Ŀ� ���۾�
'

c2_hh <- read_excel("~/Desktop/KCDC/HLAimputation/02.HLAepitope_matching/20230508_epletv2/c2_header.xlsx",sheet = 4)
c2_hh
c2_hh$V5

#c2_hh %>% t() %>% as.data.frame() -> c2_header_index 

colnames(c2_h1) <- seq(1:ncol(c2_h1))
colnames(c2_h2) <- seq(1:ncol(c2_h2))

c2_hh$V5[grepl("DR",c2_hh$V5)]
grepl("DR",c2_hh$V5)

head(c2_h1)

c2_h1 %>% select(c2_hh$V4) -> c2_h1
c2_h2 %>% select(c2_hh$V4) -> c2_h2

colnames(c2_h1) <- c2_hh$V5
colnames(c2_h2) <- c2_hh$V5


head(c2_h1)
head(c2_h2)

c2_hh$V5

c2_h1_sm <- c2_h1[,c("RecInfo","DonInfo","Nr_allDRB","Nr_AbDRB","Nr_otDRB", "Nr_allDQB","Nr_AbDQB","Nr_otDQB")]
c2_h2_sm <- c2_h1[,c("RecInfo","DonInfo","Nr_allDRB","Nr_AbDRB","Nr_otDRB", "Nr_allDQB","Nr_AbDQB","Nr_otDQB")]
head(c2_h1_sm)
head(c2_h2_sm)

c2_h1_sm %>% #merge(c2_h2_sm,by = "RecInfo") %>% #count(DonInfo.x == DonInfo.y)
  merge(c2_h2_sm,by = c("RecInfo","DonInfo")) %>% #head()
  mutate("Nr_allDRB" = ifelse(Nr_allDRB.x > Nr_allDRB.y,Nr_allDRB.x,Nr_allDRB.y)) %>% #head()
  mutate("Nr_AbDRB" = ifelse(Nr_AbDRB.x > Nr_AbDRB.y,Nr_AbDRB.x, Nr_AbDRB.y)) %>% #head()
  mutate("Nr_otDRB" = ifelse(Nr_otDRB.x > Nr_otDRB.y,Nr_otDRB.x,Nr_otDRB.y)) %>% #head()
  mutate("Nr_allDQB" = ifelse(Nr_allDQB.x > Nr_allDQB.y,Nr_allDQB.x,Nr_allDQB.y)) %>% #head()
  mutate("Nr_AbDQB" = ifelse(Nr_AbDQB.x > Nr_AbDQB.y,Nr_AbDQB.x,Nr_AbDQB.y)) %>% #head()
  mutate("Nr_otDQB" = ifelse(Nr_otDQB.x > Nr_otDQB.y,Nr_otDQB.x,Nr_otDQB.y)) %>% #head()
  select(RecInfo,DonInfo,Nr_allDRB,Nr_AbDRB,Nr_otDRB,Nr_allDQB,Nr_AbDQB,Nr_otDQB) -> c2_mm_forASO


writexl::write_xlsx(c2_mm_forASO,"SMMM/KR.KD.HLAimp.foreplet.v2_bCODE_singleMM_scoreCount.xlsx")



###### histogram
library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)
setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/hla_matchmaker")


c1 <- read_excel("c1_forWAS.xlsx")
c2 <- read_excel("c2_forWAS.xlsx")

pheno <- read_xlsx("~/Desktop/KCDC/transplantation/allogenomic/phenotype/DATASET_20250315_Final.xlsx")
head(c1)
head(c2)
head(pheno)
dim(pheno)
#pheno %>% select(SUBJNO)
pheno %>% #filter(bCODE_R %in% out$OriID) %>% 
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
    TRUE ~ as.numeric(LAST_FU_DATE - KT_DATE, units = "days"))) -> pheno

head(pheno)
head(c1)
c1 %>% left_join(pheno %>% rename(RecInfo = bCODE_R) %>% select(RecInfo,AR_TREAT,AR_BX,AR)) %>% #head
  select(RecInfo,`All Eplets`,AbVEp,OthEp,AR_TREAT,AR_BX,AR) %>% 
  pivot_longer(
    cols = c(`All Eplets`, AbVEp, OthEp),
    names_to = "Eplet_type", values_to = "Eplet_value"
  ) %>%
  pivot_longer(
    cols = c(AR_TREAT, AR_BX, AR),
    names_to = "Condition", values_to = "Status"
  ) %>%
  filter(!is.na(Eplet_value)) %>%
  mutate(
    # 보기 좋은 순서/라벨
    Condition = factor(Condition, levels = c("AR_TREAT","AR_BX","AR")),
    Status    = factor(Status, levels = c(0,1), labels = c("Control","Case"))
  ) -> c1_pheno

ggplot(c1_pheno, aes(x = Eplet_value, fill = Status)) +
  geom_histogram(position = "stack", bins = 20, color = "white") +
  facet_grid(Condition ~ Eplet_type, scales = "free") +
  scale_fill_manual(values = c("#bdbdbd", "#3182bd"), name = "Status") +
  labs(
    x = "Score", y = "Sample count",
    title = "Distributions of Eplet metrics by AR conditions (Class I)"
  ) +
  theme_bw()
  


colnames(c2)
c2 %>% mutate(All_Other_ClassII = Total_eps - All_ABV_ClassII) %>%
  left_join(pheno %>% rename(RecInfo = bCODE_R) %>% select(RecInfo,AR_TREAT,AR_BX,AR)) %>% #head
  select(RecInfo,All_ABV_ClassII:Nr_otDPA,All_Other_ClassII,AR_TREAT,AR_BX,AR) %>%  #head()
  pivot_longer(
    cols = c(All_ABV_ClassII:Nr_otDPA,All_Other_ClassII),
    names_to = "Eplet_type", values_to = "Eplet_value"
  ) %>%
  pivot_longer(
    cols = c(AR_TREAT, AR_BX, AR),
    names_to = "Condition", values_to = "Status"
  ) %>%
  filter(!is.na(Eplet_value)) %>%
  mutate(
    # 보기 좋은 순서/라벨
    Condition = factor(Condition, levels = c("AR_TREAT","AR_BX","AR")),
    Status    = factor(Status, levels = c(0,1), labels = c("Control","Case"))
  ) -> c2_pheno
head(c2_pheno)
colnames(c2)
table(c2_pheno$Eplet_type)


ggplot(c2_pheno %>% filter(Eplet_type %in% c("All_ABV_ClassII","All_Other_ClassII","Total_eps")), aes(x = Eplet_value, fill = Status)) +
  geom_histogram(position = "stack", bins = 20, color = "white") +
  facet_grid(Condition ~ Eplet_type, scales = "free") +
  scale_fill_manual(values = c("#bdbdbd", "#3182bd"), name = "Status") +
  labs(
    x = "Score", y = "Sample count",
    title = "Distributions of Eplet metrics by AR conditions (Class II)"
  ) +
  theme_bw()


## data prop

class1 <- read_xlsx("HLAmathmaker_classI_prop_20251103.xlsx")

class2 <- read_xlsx("HLAmathmaker_classII_prop_20251103.xlsx")
pirche <- read_xlsx("../PIRCHE/pirche_result_KOTRY_250421.xlsx")

head(pirche)
