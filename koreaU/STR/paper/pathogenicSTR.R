library(tidyverse)
## patho


load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.info.RData")
load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.info.RData")
trgt <- tgrt
tgrt <- NULL

head(trgt)
head(eh)


trgt %>% filter(!str_detect(TRID,"chr")) ->  trgt_guo_patho
eh %>% filter(!str_detect(STR_ID,"chr")) ->  eh_guo_patho
trgt_path <- read_table()
#write.table(trgt_guo_patho,"/Users/ksmpooh/Desktop/KU/@research/STR/patho/beforeQC/trgt_guo_patho_beforeQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(eh_guo_patho,"/Users/ksmpooh/Desktop/KU/@research/STR/patho/beforeQC/eh_guo_patho_beforeQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_simple_pass_intersect.RData")
load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_simple_pass_intersect.RData")

head(eh_simple_pass_intersect)
head(trgt_simple_pass_intersect)

trgt_simple_pass_intersect %>% filter(!str_detect(TRID,"chr")) ->  trgt_strdb_simple_patho
eh_simple_pass_intersect %>% filter(!str_detect(STR_ID,"chr")) ->  eh_strdb_simple_patho
#head(trgt_strdb_simple_patho)
#write.table(trgt_strdb_simple_patho,"/Users/ksmpooh/Desktop/KU/@research/STR/patho/beforeQC/trgt_strdb_simple_patho.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#write.table(eh_strdb_simple_patho,"/Users/ksmpooh/Desktop/KU/@research/STR/patho/beforeQC/eh_strdb_simple_patho.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#####

load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/eh_complex_pass.RData")
load("/Users/ksmpooh/Desktop/KU/@research/STR/Rdata/v2/trgt_complex_pass.RData")

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/patho/beforeQC/")
trgt_guo_patho <- read_table("trgt_guo_patho_beforeQC.txt")
eh_guo_patho <- read_table("eh_guo_patho_beforeQC.txt")
trgt_strdb_simple_patho <- read_table("trgt_strdb_simple_patho.txt")
eh_strdb_simple_patho <- read_table("eh_strdb_simple_patho.txt")



final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
head(final_ref)
final_ref %>% filter(str_detect(MOTIFS,",")) %>% select(ID,MOTIFS) %>% separate_rows(MOTIFS, sep = ",") %>%
  mutate(new_ID = paste(ID, MOTIFS, sep = "_")) -> complex_ID
colnames(complex_ID)[1] <- c("STR_ID")
eh_complex_pass %>% mutate(new_ID = ifelse(str_detect(STR_ID,"_"),STR_ID,paste(STR_ID, RU, sep = "_"))) %>%
  mutate(STR_ID = str_split_fixed(new_ID,"_",2)[,1]) -> eh_patho_merge
head(eh_patho_merge)
head(complex_ID)
head(eh_complex_pass)
head(trgt_complex_pass)
trgt_complex_pass %>% select(ID,TRID,MOTIFS,Allele,MC)
trgt_complex_pass %>% separate_rows(MOTIFS, sep = ",") %>% mutate(new_ID = paste0(TRID,"_",MOTIFS)) %>% select(-MC) -> trgt_complex_pass_MOTIFS
trgt_complex_pass %>% #separate_rows(MOTIFS, sep = ",") %>% mutate(new_ID = paste0(TRID,"_",MOTIFS)) %>% 
  separate_rows(MC, sep = "_") %>% select(MC) -> trgt_complex_pass_MC
trgt_complex_pass_MC
#trgt_complex_pass %>% separate_rows(MC, sep = ",") %>% mutate(new_ID = paste0(TRID,"_",MOTIFS)) %>% select(-MC) -> trgt_complex_pass_MOTIFS

#select(-MOTIFS) -> trgt_complex_pass_MC
head(trgt_complex_pass_MOTIFS)
head(trgt_complex_pass_MC)
dim(trgt_complex_pass_MOTIFS)
dim(trgt_complex_pass_MC)

trgt_complex_pass_MOTIFS %>% cbind(trgt_complex_pass_MC) -> trgt_path_merge



head(eh_complex_pass)
head(trgt_complex_pass)
head(trgt_guo_patho)
head(eh_guo_patho)
head(trgt_strdb_simple_patho)
head(eh_strdb_simple_patho)


### trgt
head(trgt_path_merge)
trgt_path_merge %>% select(ID,TRID,MOTIFS,Allele,MC,new_ID) %>% #head()
  pivot_wider(names_from = Allele,values_from = MC) %>% rename("MC1" = allele_1,"MC2" = allele_2) -> trgt_path_merge_gt
head(trgt_path_merge_gt)
trgt_path_merge %>% select(ID,TRID,MOTIFS,Allele,AP,new_ID) %>% #head()
  pivot_wider(names_from = Allele,values_from = AP) %>% rename("AP1" = allele_1,"AP2" = allele_2) -> trgt_path_merge_AP

trgt_path_merge %>% select(ID,TRID,MOTIFS,Allele,AM,new_ID) %>% #head()
  pivot_wider(names_from = Allele,values_from = AM) %>% rename("AM1" = allele_1,"AM2" = allele_2) -> trgt_path_merge_AM


head(trgt_path_merge_gt)
head(trgt_path_merge_AM)
head(trgt_path_merge_AP)

trgt_path_merge_gt %>% left_join(trgt_path_merge_AM) %>% left_join(trgt_path_merge_AP) -> trgt_path_merge_merge

head(trgt_path_merge_merge)
head(trgt_path_merge_merge)
trgt_path_merge_merge %>% filter(MC1 == ".")
trgt_path_merge_merge %>% filter(MC2 == ".")
trgt_path_merge_merge %>% filter(is.na(AP1))
trgt_path_merge_merge %>% filter(is.na(AP2))
trgt_path_merge_merge %>% filter(AP2 == ".")
trgt_path_merge_merge %>% filter(AP1 == ".")
trgt_path_merge_merge %>% filter(AP1 == "0")

trgt_path_merge_merge %>%
  mutate(across(AM1:AP2, ~ ifelse(. == ".", -1, .))) ->  trgt_path_merge_merge
head(trgt_path_merge_merge)

trgt_path_merge_merge %>% mutate(TRGT_STR1 = ifelse(MC1 < MC2, MC1,MC2),TRGT_STR2 = ifelse(MC1 < MC2, MC2,MC1)) %>%
  mutate(TRGT_AM1 = ifelse(MC1 < MC2, AM1,AM2),TRGT_AM2 = ifelse(MC1 < MC2, AM2,AM1)) %>%
  mutate(TRGT_AP1 = ifelse(MC1 < MC2, AP1,AP2),TRGT_AP2 = ifelse(MC1 < MC2, AP2,AP1)) %>% 
  select(-MC1,-MC2,-AM1,-AM2,-AP1,-AP2)  -> trgt_path_merge_merge_processing


trgt_path_merge_merge_processing %>% filter(ID != "NIH23J3904558") -> trgt_path_merge_merge_processing
head(trgt_path_merge_merge_processing)
colnames(trgt_path_merge_merge_processing)[2] <- "STR_ID"
trgt_path_merge_merge_processing %>% write.table("~/Desktop/KU/@research/STR/trgt/trgt_patho_divided_motif_strdb_processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")

  
  
head(eh_complex_pass)
head(eh_guo_patho)
head(eh_strdb_simple_patho)
#eh_complex_pass %>% rbind(eh_strdb_simple_patho) %>% select(-CHROM,-POS) %>% rbind(eh_guo_patho) -> eh_patho_merge


eh_patho_merge %>% select(-type,-ADSP,-ADFL,-ADIR) %>% #filter(FILTER == "PASS") %>%
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>%
  mutate(ALT = ifelse(ALT == ".",0,ALT)) %>% mutate(GT = str_split_fixed(GT,"/",2),ALT= str_split_fixed(ALT,",",2)) -> eh_patho_merge_gt

eh_patho_merge_gt %>% mutate(STR1 = ifelse(GT[,1] == "1",ALT[,1],ifelse(GT[,1] == "2",ALT[,2],REF))) %>%
  mutate(STR2 = ifelse(GT[,2] == "1",ALT[,1],ifelse(GT[,2] == "2",ALT[,2],REF))) %>% select(-ALT,-GT,-REF) %>% mutate(STR1 = as.integer(STR1),STR2 = as.integer(STR2)) -> eh_patho_merge_gt

head(eh_patho_merge)

eh_patho_merge %>% select(ID,STR_ID,type,new_ID) %>% mutate(type1 = str_split_fixed(type,"/",2)[,1],type2 = str_split_fixed(type,"/",2)[,2]) %>% 
  select(-type) -> eh_patho_merge_type

head(eh_patho_merge_gt)
head(eh_patho_merge_type)

eh_patho_merge_gt %>% select(-FILTER) %>% left_join(eh_patho_merge_type) -> eh_patho_merge_processing

eh_patho_merge_processing %>% mutate(EH_STR1 = ifelse(STR1 < STR2, STR1,STR2),EH_STR2 = ifelse(STR1 < STR2, STR2,STR1)) %>%
  mutate(EH_type1 = ifelse(STR1 < STR2, type1,type2),EH_type2 = ifelse(STR1 < STR2, type2,type1)) %>% 
  select(-STR1,-STR2,-type1,-type2)-> eh_patho_merge_processing

head(eh_patho_merge_processing)
eh_patho_merge_processing %>% filter(ID != "NIH20N2594890") -> eh_patho_merge_processing
ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)

ref %>% select(Revio,Illumina) -> ref
head(ref)
colnames(ref) <- c("ID","EH_ID")

eh_patho_merge_processing  %>% rename(EH_ID = ID) %>% left_join(ref) -> eh_patho_merge_processing


eh_patho_merge_processing %>% write.table("~/Desktop/KU/@research/STR/eh/eh_patho_divided_motif_strdb_processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")



##### patho


eh_patho_merge_processing <- read_table("~/Desktop/KU/@research/STR/eh/eh_patho_divided_motif_strdb_processing.txt")
trgt_path_merge_merge_processing <- read_table("~/Desktop/KU/@research/STR/trgt/trgt_patho_divided_motif_strdb_processing.txt")
 
head(eh_patho_merge_processing)
head(trgt_path_merge_merge_processing)

eh_patho_merge_processing %>% select()

eh_patho_merge_processing %>% select(STR_ID) %>% unique()

eh_patho_merge_processing %>% select(STR_ID,STR_ID,new_ID,ID,EH_STR1,EH_STR2) -> eh_patho_merge_processing_complexSTR_gt
trgt_path_merge_merge_processing %>% select(-MOTIFS) %>% left_join(eh_patho_merge_processing_complexSTR_gt) -> a
#write.table(a,"~/Desktop/KU/@research/STR/figure/complexSTR_prep.EH_TRGT_merged.txt",col.names = T,row.names = F,quote = F,sep = "\t")




