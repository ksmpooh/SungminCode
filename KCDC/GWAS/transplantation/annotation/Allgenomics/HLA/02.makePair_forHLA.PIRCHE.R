library(tidyverse)

setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/HLA_KHUref_AR")
df <- read.table("output_fd.txt",header = T)
df <- read.table("output_td.txt",header = T)

ref <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/00.sampleInfo/ALL(2019and2020).sampleID.withbCODE.xlsx",)
#ref_sev < -read.table("~/Des")
head(ref)
qcin_bCODE <- readxl::read_xlsx("~/Desktop/KCDC/transplantation/00.sampleInfo/KOTRY.open.QCed.list.xlsx")
head(qcin_bCODE)
ref %>% filter(bCODE %in% qcin_bCODE$bCODE) -> ref
head(ref)
head(ref)
#A*	A*	B*	B*	C*	C*	DRB1*	DRB1*	DRB3*	DRB3*	DRB4*	DRB4*	DRB5*	DRB5*	DQA1*	DQA1*	DQB1*	DQB1*	DPA1*	DPA1*	DPB1*	DPB1*
df %>% mutate(ID = str_split_fixed(ID,"=",2)[,1]) -> df
head(df)
head(ref)

ref %>% filter(KBA_ID %in% df$ID) %>% filter(type == "KR") %>% 
  merge(ref %>% filter(KBA_ID %in% df$ID) %>% filter(type == "KD"),all = T,by='ref') %>% na.omit() -> pair_ref
head(ref)
head(df)
head(pair_ref)
pair_ref %>% select(KBA_ID.x,KBA_ID.y) %>% pivot_longer(KBA_ID.x:KBA_ID.y) %>% select(value) %>% rename(ID = value) %>%
  left_join(ref %>% rename(ID = KBA_ID) %>% select(ID,OriID)) %>% left_join(df) -> a
head(a)

colnames(a)

a[,c("OriID","HLA_A.1","HLA_A.2","HLA_B.1","HLA_B.2","HLA_C.1","HLA_C.2","HLA_DRB1.1","HLA_DRB1.2","HLA_DQA1.1","HLA_DQA1.2","HLA_DQB1.1","HLA_DQB1.2","HLA_DPA1.1","HLA_DPA1.2","HLA_DPB1.1","HLA_DPB1.2")] -> a
head(a)

a$empty <- ""
"empty"
#A*	A*	B*	B*	C*	C*	DRB1*	DRB1*	DRB3*	DRB3*	DRB4*	DRB4*	DRB5*	DRB5*	DQA1*	DQA1*	DQB1*	DQB1*	DPA1*	DPA1*	DPB1*	DPB1*
a[,c("OriID","HLA_A.1","HLA_A.2","HLA_B.1","HLA_B.2","HLA_C.1","HLA_C.2","HLA_DRB1.1","HLA_DRB1.2","empty","empty","empty","empty","empty","empty","HLA_DQA1.1","HLA_DQA1.2","HLA_DQB1.1","HLA_DQB1.2","HLA_DPA1.1","HLA_DPA1.2","HLA_DPB1.1","HLA_DPB1.2")] -> out

colnames_str <- "Patient_Donor_ID A* A* B* B* C* C* DRB1* DRB1* DRB3* DRB3* DRB4* DRB4* DRB5* DRB5* DQA1* DQA1* DQB1* DQB1* DPA1* DPA1* DPB1* DPB1*"
head(out)
# 문자열을 공백을 기준으로 분할하여 벡터로 변환
colnames_vector <- strsplit(colnames_str, " ")[[1]]
colnames(out) <- colnames_vector
head(out)


nih <- out


#write.csv(out,"~/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/KOTRY.AR.HLAimp_forPIRCHE.csv",row.names = F,quote = F)
write.csv(out,"~/Desktop/KCDC/transplantation/allogenomic/2025_AR/HLA/KOTRY.AR.HLAimp_forPIRCHE.td.csv",row.names = F,quote = F)

df %>% filter(!str_detect(ID,"NIH")) -> df_ys
ref_ys_pair <- read.table("~/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/YSprod_pair_ID.txt",header = T)

head(df_ys)
head(ref_ys_pair)
ref_ys_pair %>%

