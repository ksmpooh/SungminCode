setwd("~/Desktop/KCDC/transplantation/")
ref <- readxl::read_xlsx("00.sampleInfo/ALL(2019and2020).sampleID.withbCODE.xlsx")
head(ref)


ys_check <- c("KR01694","KR01138","KR01107","KR03892","KR04394","KR04393","KR00514","KR08118","KR04210","KR00176")

ref %>% filter(OriID %in% ys_check) -> a
#write.table(a,"allogenomic/2025_AR/check.txt",col.names = F,row.names = F,quote = F)

ys_list <- readxl::read_xlsx("allogenomic/2025_AR/sample_info/2_Main data_250122_NoDSA.xlsx",sheet = 3)
full_ys_list <- readxl::read_xlsx("allogenomic/2025_AR/sample_info/2_Main data_250122_NoDSA.xlsx",sheet = 2)
head(ys_list)
head(full_ys_list)
tail(full_ys_list)

full_ys_list %>% filter(TOTAL_ID %in% ys_list$TOTAL_ID)
full_ys_list %>% filter(!(TOTAL_ID %in% ys_list$TOTAL_ID)) %>% filter(!is.na(KOTRY_ID)) %>%
  head()

13+18
tail(full_ys_list)
head(ref)
full_ys_list %>% filter(!is.na(KOTRY_ID) | !is.na(SEV_ID)) %>% 
  mutate(prod_by = ifelse(TOTAL_ID %in% ys_list$TOTAL_ID,"YS","NIH")) -> new_df

new_df %>% filter(prod_by == "NIH") -> nih
new_df %>% filter(prod_by != "NIH") -> ys

nih %>%
  left_join(ref %>% select(KBA_ID,ref,OriID) %>% rename(KOTRY_ID = OriID,KR = KBA_ID)) %>% #head()
  left_join(ref %>% filter(type == "KD") %>% select(KBA_ID,ref) %>% rename(KD = KBA_ID)) %>% select(-ref)-> nih

head(nih)  
head(ys)
head(ys_list)
colnames(ys_list)
ys %>% left_join(ys_list) %>% select(!`...6`) -> ys
colnames(ys) <- colnames(nih)
head(ys)  
head(nih)

cel_list <- read.table("allogenomic/2025_AR/KBA_QC/all.cellID.KRKD.txt")
head(cel_list)
cel_list %>% mutate(ID = ifelse(grepl("NIH",V1),str_split_fixed(str_split_fixed(V1,"_",6)[,6],"\\.",2)[,1],0)) %>% filter(ID == 0)
cel_list %>% mutate(ID = str_split_fixed(str_split_fixed(V1,"_",6)[,6],"\\.",2)[,1]) %>% rename(cel_file_name = V1) -> cel_list

tail(cel_list)
tail(ys)
tail(cel_list)
cel_list %>% filter(ID == "0811-80-1")
ys %>% mutate(KR = str_replace_all(KR,"_","-")) %>% mutate(KD = str_replace_all(KD,"_","-")) %>% rbind(nih) -> df
head(df)

df %>% left_join(cel_list %>% rename(cel_file_name.KR = cel_file_name,KR = ID)) %>%
  left_join(cel_list %>% rename(cel_file_name.KD = cel_file_name,KD = ID)) %>%
  write.table("allogenomic/2025_AR/KBA_QC/sampleInfo.requ.byYS.20250131.txt",col.names = T,row.names = F,quote = F,sep = "\t")


