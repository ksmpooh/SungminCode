library(tidyverse)
library(stringr)
setwd("~/Desktop/KCDC/HLAimputation/MakeReferencePanel/")


df <-read.table("HLAtype_check/nomane_log/HLA.type.result.8genes.merged_529sample_forMAKEreference_ori.txt",header = T)
ref <- readxl::read_xls("~/Desktop/KCDC/????????????????????????????????????/??????????????????????????????/HLAtyping.265pairtable_with(QC.typing)_20211028.xls")
head(ref)
colnames(ref)
ref %>% select(KBA_ID.2019,KBA_ID.2020,HLAID.2019,HLAID.2020) %>% dim()

ref %>% select(KBA_ID.2019,HLAID.2019) %>% dim()
ref %>% select(KBA_ID.2020,HLAID.2019) %>% dim()


ref1 <- ref %>% select(KBA_ID.2019,HLAID.2019)
ref2 <- ref %>% select(KBA_ID.2020,HLAID.2020)
ref %>% select(KBA_ID.2019,HLAID.2019) %>% dim()
ref %>% select(KBA_ID.2020,HLAID.2019) %>% dim()

colnames(ref1) <- c("ID","HLAID")
colnames(ref2) <- c("ID","HLAID")
ref <-rbind(ref1,ref2)



#df <-read.table("HLA.4digit.520sample.nomenclean_withheader.chped",header = T)
#df <-read.table("")
head(df)

nomen <- read.table("HLAtype_check/nomane_log/HLA.4digit.520sample.nomenclean.chped")
head(nomen)
colnames(nomen) <- colnames(df)


df %>% filter(IID %in% nomen$IID) %>% select(-IID,-pID,-mID,-SEX,-PHENO) %>% 
  pivot_longer(2:17,names_to = 'gene',values_to = 'NGS.type') -> df1
head(df1)

nomen %>%  select(-IID,-pID,-mID,-SEX,-PHENO) %>% 
  pivot_longer(2:17,names_to = 'gene',values_to = 'Nomen.type') %>% #mutate("theme" = "nomen") %>%
  inner_join(df1) %>% mutate(NGS.type = paste0(str_split_fixed(NGS.type,":",3)[,1],":",str_split_fixed(NGS.type,":",3)[,2])) %>%
  mutate(gene = str_replace_all(gene,"NGS_","")) %>%
  mutate(gene = str_replace_all(gene,"\\.1","")) %>%
  mutate(gene = str_replace_all(gene,"\\.2","")) %>%
  filter(NGS.type != ":") -> df2

head(df2)  


## check nomen vs ngs

head(df2)
head(ref)
ref %>% mutate(HLAID = str_replace_all(HLAID,"H","CDC")) -> ref
colnames(df2)[1] <- "ID"
df2 %>% mutate("check" = ifelse(Nomen.type == NGS.type,"yes","no")) %>% #head()#count(check)
  filter(check == "no") %>% left_join(ref) %>% select(-ID) %>%
  select(HLAID,gene,NGS.type,check) %>% 
  arrange(HLAID,gene) %>% 
  writexl::write_xlsx("HLAtype_check/HLAtype.check.4digit.withID.xlsx")
   #%>% mutate(gene = str_replace_all(gene,"NGS_","")) %>%
  
head(df2)

df2 %>% count(gene,NGS.type) %>% #head()
  group_by(gene) %>% 
  mutate(freq = prop.table(n)) %>% inner_join(df2 %>% select(-FID) %>% unique()) %>% #head()
  head()
  #mutate("check" = ifelse(Nomen.type == NGS.type,"yes","no")) %>%
  arrange(gene,-freq) %>%
  writexl::write_xlsx("HLAtype_check/HLAtype.compare.HLA_TAPAS.xlsx")
  
  #mutate("check" = ifelse(Nomen.type == NGS.type,"yes","no")) %>% count(check)


#### 8 digit
head(df)
#nomen <-read.table("HLAtype_check/")
imgt3320 <-  read.table("HLAtype_check/HLA_ALLELE_TABLE.imgt3320.hat",header = T)
head(imgt3320)

imgt3320 %>% mutate(Nomen.type=paste0(HLA,"*",STANDARD)) %>%
  select(AlleleID,Nomen.type) -> imgt3320_v1


imgt3320_v1 %>% filter(Nomen.type %in% df1$NGS.type) %>% dim()
colnames(imgt3320_v1) <- c("IPDID_inIMGT3320","NGS.type")

df %>%
  filter(IID %in% nomen$IID) %>% select(-IID,-pID,-mID,-SEX,-PHENO) %>% 
  pivot_longer(2:17,names_to = 'gene',values_to = 'NGS.type') %>% 
  mutate(gene = str_replace_all(gene,"NGS_","")) %>%
  mutate(gene = str_replace_all(gene,"\\.1","")) %>%
  mutate(gene = str_replace_all(gene,"\\.2","")) %>%
  #inner_join(imgt3320_v1) %>% #count(is.na(IPD_ID))
  count(gene,NGS.type) %>% 
  mutate(freq = prop.table(n)) %>%
  left_join(imgt3320_v1) %>% select(-n) %>% #head()
  mutate(IPDID_inIMGT3320 = ifelse(is.na(IPDID_inIMGT3320),"NA",IPDID_inIMGT3320)) %>%
  arrange(gene,-freq) %>% 
  writexl::write_xlsx("HLAtype_check/HLAtype.compare.HLA_TAPAS_IMGT3320.8digit.xlsx")
  



head(imgt3320_v1)
df2 %>% count(gene,NGS.type) %>% #head()
  group_by(gene) %>% 
  mutate(freq = prop.table(n)) %>% #head()
  left_join(imgt3320_v1) %>% select(-n) %>% count(IPDID_inIMGT3320) #,head()
  inner_join(df2 %>% select(-FID) %>% unique()) %>% #head()
  head()
#mutate("check" = ifelse(Nomen.type == NGS.type,"yes","no")) %>%
arrange(gene,-freq) %>%
  writexl::write_xlsx("HLAtype_check/HLAtype.compare.HLA_TAPAS.xlsx")



## with ID
head(imgt3320_v1)
head(imgt3320)


df %>%
  filter(IID %in% nomen$IID) %>% select(-IID,-pID,-mID,-SEX,-PHENO) %>% 
  pivot_longer(2:17,names_to = 'gene',values_to = 'NGS.type') %>% 
  mutate(gene = str_replace_all(gene,"NGS_","")) %>%
  mutate(gene = str_replace_all(gene,"\\.1","")) %>%
  mutate(gene = str_replace_all(gene,"\\.2","")) %>% #head()
  left_join(imgt3320_v1) %>% #head()
  filter(is.na(IPDID_inIMGT3320)) %>% mutate(ID = FID) %>%
  mutate(IPDID_inIMGT3320 = ifelse(is.na(IPDID_inIMGT3320),"NA",IPDID_inIMGT3320)) %>%
  left_join(ref) %>% #head()
  select(HLAID,gene,NGS.type,IPDID_inIMGT3320) %>% 
  arrange(HLAID,gene) %>%
  writexl::write_xlsx("HLAtype_check/HLAtype.check.8digit.withID.xlsx")


