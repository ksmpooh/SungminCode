### main figure 4 data processing
#mapq bamQ APscore depth 
library(tidyverse)
library(ggplot2)

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro
final_ref_pro %>% filter(!str_detect(MOTIFS,",")) %>% filter(STR_ID %in%common_STR$STR_DB) -> simple_STR

final_ref %>% filter(str_detect(MOTIFS,",")) %>% 
  separate_rows(MOTIFS, sep = ",") %>%
  group_by(chrom, start, end, ID) %>%
  mutate(priority = row_number()) %>%
  ungroup() %>%
  mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX"))) %>%
  arrange(chrom, start) %>%
  rename(STR_ID = ID) %>%
  mutate(new_ID = paste0(STR_ID,"_",MOTIFS)) -> final_ref_complex

final_ref %>% filter(str_detect(MOTIFS,"N")) -> final_ref_patho_withN
#final_ref %>% filter(str_detect(MOTIFS,"R"))



head(final_ref_pro)
head(final_ref)
final_ref %>% filter(chrom == "chrX") -> final_ref_chrX
head(final_ref_chrX)
head(concordance_bySTR_simpleSTR.afterchrXQC)
head(concordance_bySTR_simpleSTR.afterchrXQC)

concordance_bySTR_simpleSTR.afterchrXQC %>% filter(STR_ID == "RFC1")
concordance_bySTR_simpleSTR.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.concordance_bySTR_simpleSTR.afterchrXQC.txt")
concordance_bySTR_simpleSTR.afterchrXQC %>% filter(STR_ID == "RFC1")
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1")
head(concordance_bySTR_simpleSTR.afterchrXQC)
concordance_bySTR_simpleSTR.afterchrXQC %>% filter(concordance_rate == 0) -> str_concordance0
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(TRGT_STR == 0) %>% count(STR_ID) -> trgt_str0
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% 
  filter(STR_ID %in% str_concordance0$STR_ID) %>% group_by(STR_ID) %>% group_by(STR_ID) %>% summarise(meanAP = mean(TRGT_AP)) -> str_concordance0_apmean

head(str_concordance0_apmean)
str_concordance0_apmean %>% filter(meanAP < 0.7)
head(trgt_str0)
trgt_str0 %>% arrange(-n)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)

trgt_cov <- read_table("~/Desktop/KU/@research/STR/trgt/oribam_cover/trgt_normal_patho_oribam_coverage.txt") %>% select(-STR_region,-type,-platform) %>% 
  filter(STR_ID %in% simple_STR$STR_ID)
head(trgt_cov)
#trgt_cov %>% filter(STR_ID == "RFC1")
trgt_cov %>% select(ID) %>% unique() %>% dim()

#colnames(trgt_cov) <- c("TRGT_meandepth","TRGT_meanbaseq","TRGT_meanmapq","ID","STR_ID") 
colnames(trgt_cov) <- c("meandepth","meanbaseq","meanmapq","ID","STR_ID") 

eh_cov <- read_table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_oribam_coverage.txt") %>% select(meandepth:ID) %>% filter(STR_ID %in% simple_STR$STR_ID)
#eh_cov_patho <- read_table("~/Desktop/KU/@research/STR/eh/") #%>% select(meandepth:ID) %>% filter(STR_ID %in% common_STR$STR_DB)

eh_cov %>% select(ID) %>% unique() %>% dim()
head(eh_cov)
#colnames(eh_cov) <- c("eh_meandepth","eh_meanbaseq","eh_meanmapq","STR_ID","ID") 
colnames(eh_cov) <- c("meandepth","meanbaseq","meanmapq","STR_ID","ID") 


head(eh_cov)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1")
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(TRGT_STR)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(TRGT_AP == 0) %>% filter(STR_ID != "RFC1") %>%
  left_join(trgt_cov %>% rename(trgt_meandepth = meandepth,trgt_meanbaseq = meanbaseq,trgt_meanmapq = meanmapq)) %>%
  left_join(eh_cov) %>% rename(eh_meandepth = meandepth,eh_meanbaseq = meanbaseq,eh_meanmapq = meanmapq) %>%
  write.table("~/Desktop/KU/@research/STR/figure/discussion/eh_trgt_merge_simple_pass_intersect_forConcordance.AP_0.INFO.afterchrXQC",col.names = T,row.names = F,quote = F,sep = "\t")

head(trgt_cov)
head(eh_cov)

head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC 
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>% group_by(ID) %>% summarise(AP_mean = mean(TRGT_AP)) -> AP_mean_byID_simpleSTR
head(AP_mean_byID_simpleSTR)
write.table(AP_mean_byID_simpleSTR,"~/Desktop/KU/@research/STR/figure/figure4/f4.AP_mean_byID_simpleSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)


eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID) %>% filter(STR_ID != "RFC1") %>%
  summarise(concordance_rate = mean(check)) %>% write.table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/concordance_rate_byID.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
   
  

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% count(AP0.5) %>% mutate(prop=prop.table(n))
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID,AP0.5) %>% summarise(concordance_rate = mean(check)) -> concordacne_byID_AP_simpleSTR_afterchrXQC
head(concordacne_byID_AP_simpleSTR_afterchrXQC)
#write.table(concordacne_byID_AP_simpleSTR_afterchrXQC,"~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_simpleSTR_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    TRUE ~ '[50~Inf)')) %>% #count(STR_length,AP0.5) %>% mutate(prop=prop.table(n))
  group_by(ID,STR_length,AP0.5) %>% summarise(concordance_rate = mean(check)) -> concordacne_byID_AP0.5_str_length_simpleSTR_afterchrXQC
#head(concordacne_byID_AP_simpleSTR_afterchrXQC)
write.table(concordacne_byID_AP0.5_str_length_simpleSTR_afterchrXQC,"~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_str_length_simpleSTR_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")



eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% 
  mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% #head()
  group_by(STR_ID,AP0.5) %>% summarise(concordance_rate = mean(check)) -> concordacne_bySTR_AP_simpleSTR_afterchrXQC

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% 
  left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~100)",
    TRUE ~ '[100~Inf)')) %>% 
  mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% #count(STR_length,AP0.5) %>% mutate(prop=prop.table(n))
  group_by(ID,STR_length,AP0.5) %>%  
  summarise(concordance_rate = mean(check)) -> concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC
write.table(concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC,"~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

head(concordacne_byID_AP_simpleSTR_afterchrXQC)
head(concordacne_bySTR_AP_simpleSTR_afterchrXQC)




concordacne_byID_AP_simpleSTR_afterchrXQC
concordacne_byID_AP_simpleSTR_afterchrXQC %>% ggplot(aes(x=AP0.5,y=concordance_rate)) + 
  geom_violin()

concordacne_bySTR_AP_simpleSTR_afterchrXQC
concordacne_bySTR_AP_simpleSTR_afterchrXQC %>% ggplot(aes(x=AP0.5,y=concordance_rate)) +
  geom_hex(bins=15)


eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% 
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% 
  #mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% #head()
  group_by(STR_ID,TRGT_AP) %>% summarise(concordance_rate = mean(check)) %>%
  ggplot(aes(x=TRGT_AP,y=concordance_rate)) +
  geom_hex(bins=15)


  
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  select(ID:EH_STR,TRGT_AP) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% #left_join(final_ref_pro) %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% count(TRGT_AP) -> a
  
a %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.STR.allele.count.byAP.txt",col.names = T,row.names = F,quote = F,sep = "\t")

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  select(ID:EH_STR,TRGT_AP) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% #left_join(final_ref_pro) %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% count(TRGT_AP,check) -> b
  
b %>% group_by(TRGT_AP) %>%
  mutate(prop = prop.table(n)) %>%
  write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.STR.allele.count.match.byAP.txt",col.names = T,row.names = F,quote = F,sep = "\t")


eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  select(ID:EH_STR,TRGT_AP) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% #left_join(final_ref_pro) %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% group_by(ID,TRGT_AP) %>% 
  summarise(concordance_rate = mean(check)) -> c

c %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.sample.concordance.rate.byAP.txt",col.names = T,row.names = F,quote = F,sep = "\t")



eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  select(ID:EH_STR,TRGT_AP) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>%
  group_by(ID,STR_length,TRGT_AP) %>% 
  summarise(concordance_rate = mean(check)) -> d

d %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.sample.concordance.rate.bySTR_length_AP.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  select(ID:EH_STR,TRGT_AP) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>%
  mutate(STR_length = TRGT_STR * RU.length) %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>%
  group_by(STR_length,TRGT_AP) %>% 
  count(check) -> e
head(e)
e %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.match.count.bySTR_length_AP.txt",col.names = T,row.names = F,quote = F,sep = "\t")

## ap by simple patho simple_patho complex
## RFC1 -> complex
patho_complex_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")
head(patho_complex_processing)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1 <- eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1")
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
head(final_ref_patho_withN)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID %in% final_ref_patho_withN$ID) -> path_withN
head(path_withN)
path_withN %>% count(TRGT_STR == EH_STR)
path_withN %>% filter(STR_ID == "ARX_2") %>% count(TRGT_STR,EH_STR)
path_withN %>% count(STR_ID,TRGT_STR,EH_STR)
 
head(patho_complex_processing)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1 <- eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1 %>% mutate(new_ID = STR_ID,STR_type = "complex")
#patho_complex_processing %>% group_by(STR_ID) %>% summarise(TRGT_AP = mean(AP)) -> a
  

patho_complex_processing %>% rbind(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1 %>% rename(AP = TRGT_AP,AM = TRGT_AM)) %>%
  group_by(ID,STR_ID,allele) %>% summarise(TRGT_AP = mean(AP)) %>% 
  mutate(type = "complexSTR",patho = "Pathogenic") -> ap_complexSTR
patho_complex_processing 

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(ID,STR_ID,TRGT_AP,allele) %>% filter(STR_ID != "RFC1") %>%
  mutate(TRGT_AP = ifelse(TRGT_AP == -1,0,TRGT_AP)) %>% 
  mutate(type = "simpleSTR",patho = ifelse(str_detect(STR_ID,"chr"),"Normal","Pathogenic")) %>%  
  #filter(patho == "Pathogenic") %>% group_by(STR_ID) %>%summarise(TRGT_AP = mean(TRGT_AP)) -> a
  rbind(ap_complexSTR) %>% group_by(ID) %>% summarise(APmean = mean(TRGT_AP)) -> APmean_byID_allSTR

head(APmean_byID_allSTR)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1")
#write.table(APmean_byID_allSTR,"~/Desktop/KU/@research/STR/figure/figure4/f4.APmean_byID_allSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1") %>% count(TRGT_STR==EH_STR)
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(str_detect(MOTIFS)) %>% count(TRGT_STR==EH_STR)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(ID,STR_ID,TRGT_AP,allele) %>% filter(STR_ID != "RFC1") %>%
  mutate(TRGT_AP = ifelse(TRGT_AP == -1,0,TRGT_AP)) %>%
  mutate(type = "simpleSTR",patho = ifelse(str_detect(STR_ID,"chr"),"Normal","Pathogenic")) %>% 
  rbind(ap_complexSTR) %>% group_by(ID,type,patho) %>% summarise(APmean = mean(TRGT_AP)) -> APmean_byID_bySTRtype
head(APmean_byID_bySTRtype)
#write.table(APmean_byID_bySTRtype,"~/Desktop/KU/@research/STR/figure/figure4/f4.APmean_byID_bySTRtype.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

APmean_byID_allSTR
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
#APmean_byID_bySTRtype %>% 
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>%
  mutate(TRGT_AP = ifelse(TRGT_AP == -1,0,TRGT_AP)) %>%
  left_join(final_ref_pro) %>% 
  mutate(STR_length = TRGT_STR*RU.length,check = ifelse(TRGT_STR == EH_STR,1,0)) %>% 
  select(ID,STR_ID,STR_length,TRGT_AP,check,allele)

head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1)

patho_complex_processing %>% left_join(final_ref_complex %>% mutate(RU.length = str_length(MOTIFS)) %>% select(new_ID,RU.length)) %>%
  rbind(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1 %>% 
          rename(AP = TRGT_AP,AM = TRGT_AM) %>% left_join(final_ref %>% mutate(RU.length = str_length(MOTIFS)) %>% select(ID,RU.length) %>% rename(STR_ID = ID))) %>%
  mutate(STR_length = TRGT_STR*RU.length,check = ifelse(TRGT_STR == EH_STR,1,0)) %>% #head()
  group_by(ID,STR_ID,allele) %>%
  reframe(STR_length = sum(STR_length),TRGT_AP = mean(AP),check = mean(check)) %>% 
  mutate(check = ifelse(check == 1,1,0)) %>%
  rbind(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
          mutate(TRGT_AP = ifelse(TRGT_AP == -1,0,TRGT_AP)) %>%
          left_join(final_ref_pro) %>% 
          mutate(STR_length = TRGT_STR*RU.length,check = ifelse(TRGT_STR == EH_STR,1,0)) %>% 
          select(ID,STR_ID,STR_length,TRGT_AP,check,allele)) -> AP_STRlength_all_STRallele

head(AP_STRlength_all_STRallele)

AP_STRlength_all_STRallele
AP_STRlength_all_STRallele %>% filter(STR_ID =="RFC1")
AP_STRlength_all_STRallele %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% count(TRGT_AP) -> a

a %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.STR.allele.count.byAP.txt",col.names = T,row.names = F,quote = F,sep = "\t")

AP_STRlength_all_STRallele %>% 
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% count(TRGT_AP,check) -> b

b %>% group_by(TRGT_AP) %>%
  mutate(prop = prop.table(n)) %>%
  write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.STR.allele.count.match.byAP.txt",col.names = T,row.names = F,quote = F,sep = "\t")


AP_STRlength_all_STRallele %>% 
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% group_by(ID,TRGT_AP) %>% 
  summarise(concordance_rate = mean(check)) -> c

c %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.sample.concordance.rate.byAP.txt",col.names = T,row.names = F,quote = F,sep = "\t")



AP_STRlength_all_STRallele %>%
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>%
  group_by(ID,STR_length,TRGT_AP) %>% 
  summarise(concordance_rate = mean(check)) -> d

d %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.sample.concordance.rate.bySTR_length_AP.txt",col.names = T,row.names = F,quote = F,sep = "\t")


AP_STRlength_all_STRallele %>% 
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    STR_length < 100 ~ "[50~100)",
    STR_length < 150 ~ "[100~150)",
    TRUE ~ '[150~Inf)')) %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>%
  group_by(STR_length,TRGT_AP) %>% 
  count(check) -> e
head(e)
e %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/f4.match.count.bySTR_length_AP.txt",col.names = T,row.names = F,quote = F,sep = "\t")

  
head(final_ref_complex)
head(AP_STRlength_all_STRallele)

AP_STRlength_all_STRallele %>% mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% count(AP0.5) %>% mutate(prop=prop.table(n))
AP_STRlength_all_STRallele %>% mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% 
  group_by(ID,AP0.5) %>% summarise(concordance_rate = mean(check)) -> concordacne_byID_AP_allSTR_afterchrXQC
#head(concordacne_byID_AP_simpleSTR_afterchrXQC)
write.table(concordacne_byID_AP_allSTR_afterchrXQC,"~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_allSTR_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

AP_STRlength_all_STRallele %>% mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% 
  mutate(STR_length = case_when(
    STR_length < 50 ~ "[0~50)",
    TRUE ~ '[50~Inf)')) %>% #count(STR_length,AP0.5) %>% mutate(prop=prop.table(n))
  group_by(ID,STR_length,AP0.5) %>% summarise(concordance_rate = mean(check)) -> concordacne_byID_AP0.5_str_length_allSTR_afterchrXQC
#head(concordacne_byID_AP0.5_str_length_allSTR_afterchrXQC)
write.table(concordacne_byID_AP0.5_str_length_allSTR_afterchrXQC,"~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_str_length_allSTR_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")

head(concordacne_bySTR_AP_allSTR_afterchrXQC)

AP_STRlength_all_STRallele %>% 
  mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% #head()
  group_by(STR_ID,AP0.5) %>% summarise(concordance_rate = mean(check)) -> concordacne_bySTR_AP_allSTR_afterchrXQC

AP_STRlength_all_STRallele %>% 
  mutate(STR_length = case_when(
    STR_length < 100 ~ "[0~100)",
    TRUE ~ '[100~Inf)')) %>% 
  mutate(AP0.5 = ifelse(TRGT_AP >0.5,"(0.5~1]","[0~0.5]")) %>% #count(STR_length,AP0.5) %>% mutate(prop=prop.table(n))
  group_by(ID,STR_length,AP0.5) %>%  
  summarise(concordance_rate = mean(check)) -> concordacne_byID_AP0.5_strlength_allSTR_afterchrXQC
write.table(concordacne_byID_AP0.5_strlength_simpleSTR_afterchrXQC,"~/Desktop/KU/@research/STR/figure/figure4/f4.concordacne_byID_AP0.5_strlength_allSTR_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")



#################### bam QC 다시 확인해야함!!!!!!!!!!!!!
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/sup.figure/concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC")

concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC_RFC1 <- concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% filter(STR_ID == "RFC1")
head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)

#concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>%


#merge_coverage <- read_table("~/Desktop/KU/@research/STR/figure/extra_info/STR.mean.bam.stats.txt")
head(merge_coverage)
head(df)
dim(df)
head(merge_coverage)
head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% 
  mutate(concordance_rate_range = case_when(
    concordance_rate < 0.4 ~ "[0, 0.4)",
    concordance_rate < 0.8 ~ "[0.4, 0.8)",
    TRUE ~ "[0.8, 1]"
  )) %>%
  ggplot(aes(x=trgt_meanmapq,y=eh_meanmapq)) +
  geom_point(alpha=0.5) + 
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0,60)) + 
  labs(x="mean mapQ (LRS)",y="mean mapQ (SRS)") + 
  facet_grid(~factor(concordance_rate_range,levels= c("[0, 0.4)","[0.4, 0.8)","[0.8, 1]")))

# 굿
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>%
  pivot_longer(trgt_meanmapq:eh_meanmapq) %>% #head()
  ggplot(aes(x=value,fill=name)) + 
  geom_density()

concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% select(STR_ID,concordance_rate,trgt_meanmapq,eh_meanmapq) %>%
  filter(concordance_rate == 0) %>%
  #pivot_longer(trgt_meanmapq:eh_meanmapq) %>% #head()
  ggplot(aes(x=trgt_meanmapq,y=eh_meanmapq)) + 
  geom_point()


head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% 
  select(STR_ID,concordance_rate,trgt_meanbaseq,eh_meanbaseq) %>%
 # filter(concordance_rate == 0) %>%
  pivot_longer(trgt_meanbaseq:eh_meanbaseq) %>% #head()
  ggplot(aes(x=value,fill=name)) + 
  geom_density()
head(concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC)
concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% select(STR_ID,concordance_rate,trgt_meandepth,eh_meandepth) %>%
  filter(concordance_rate == 0) %>%
  pivot_longer(trgt_meandepth:eh_meandepth) %>% #head()
  ggplot(aes(x=value,fill=name)) + 
  geom_density()


concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC %>% #ㅇhead()
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "0",
    concordance_rate < 0.5 ~ "(0, 0.5)",
    concordance_rate < 1 ~ "[0.5, 1)",
    TRUE ~ "1"
  )) %>% #head()
  select(-trgt_meandepth,-eh_meandepth) %>% 
  pivot_longer(trgt_meanbaseq:eh_meanmapq) %>% #head()
  mutate(platform = ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  mutate(g = ifelse(str_detect(name,"meanbaseq"),'baseQ','mapQ')) %>% #head()
  ggplot(aes(x=value,fill=platform)) + 
  geom_density() + 
  facet_grid(g~factor(concordance_rate_range,levels= c("0","(0, 0.5)","[0.5, 1)","1")),scales = "free")
  
concordance_bySTR_simpleSTR.afterchrXQC %>% pivot_longer(cols = trgt_meandepth:eh_meanmapq) %>%  #head()
  mutate(g = ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%  #head()
  group_by(name,g,concordance_rate_range) %>%
  summarise(mean = mean(value)) %>%  head()
  
######## ektl
patho_complex_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt") %>% 
  filter(STR_ID != "RFC1")

normalSTR.STRlength100.with.bamstats <- read_table("~/Desktop/KU/@research/STR/figure/figure4/normalSTR.STRlength100.with.bamstats.txt")

pathoSTR.STRlength.with.bamstats <- read_table("~/Desktop/KU/@research/STR/figure/figure4/pathoSTR.STRlength.with.bamstats.txt")
head(pathoSTR.STRlength.with.bamstats)
head(normalSTR.STRlength100.with.bamstats)
eh_LC_normal <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/eh/eh_simpleSTR_LC.txt")
head(eh_LC_normal)

pathoSTR.STRlength.with.bamstats %>% filter(STR_ID != "RFC1") %>% #head()
  filter(TRGT_STR_length >= 100 & EH_STR_length >= 100) %>% 
  rbind(normalSTR.STRlength100.with.bamstats %>% mutate(check = ifelse(TRGT_STR == EH_STR,1,0),type = 'normal') %>%
          select(-TRGT_STR,-EH_STR,-RU.length)) -> a


#a %>%
head(a)
a %>% 
  ggplot(aes(x=trgt_meanbaseq,y= eh_meanbaseq,color=check)) + 
  geom_point()

a %>% 
  ggplot(aes(x=trgt_meanmapq,eh_meanmapq,color=check)) + 
  geom_point()


a %>% mutate(length = ifelse(TRGT_STR_length >= EH_STR_length,"trgt","eh")) %>%
  ggplot(aes(x=trgt_meandepth,eh_meandepth,color=factor(check))) + 
  geom_point(alpha = 0.6) + 
  coord_equal() + 
  facet_grid(~length)

a %>% mutate(length = ifelse(TRGT_STR_length >= EH_STR_length,"trgt","eh")) %>%
  ggplot(aes(x=trgt_meanmapq,eh_meanmapq,color=factor(check))) + 
  geom_point(alpha = 0.6) + 
  coord_equal() + 
  facet_grid(~length)

a %>% mutate(length = ifelse(TRGT_STR_length >= EH_STR_length,"trgt","eh")) %>%
  ggplot(aes(x=trgt_meanbaseq,eh_meanbaseq,color=factor(check))) + 
  geom_point(alpha = 0.6) + 
  coord_equal() + 
  facet_grid(~length)

a %>% left_join(eh_LC_normal %>% select(STR_ID,ID,LC)) %>% filter(TRGT_STR_length < EH_STR_length) %>% count(trgt_meandepth>=eh_meandepth)
'
`trgt_meandepth >= eh_meandepth`     n
<lgl>                            <int>
  1 FALSE                             4402
2 TRUE                              1293
'

a %>% left_join(eh_LC_normal %>% select(STR_ID,ID,LC)) %>% filter(TRGT_STR_length < EH_STR_length) %>% count(trgt_meandepth>=LC)

a %>% left_join(eh_LC_normal %>% select(STR_ID,ID,LC)) %>% filter(TRGT_STR_length < EH_STR_length) %>% count(trgt_meandepth>=eh_meandepth)
  #mutate(length = ifelse(TRGT_STR_length >= EH_STR_length,"trgt","eh")) %>% count(length) 
  ggplot(aes(x=trgt_meandepth,LC,color=factor(check))) + 
  geom_point(alpha = 0.6)
  #coord_equal() + 
  #facet_grid(~length)
'
  `trgt_meandepth >= LC`     n
  <lgl>                  <int>
1 FALSE                   5319
2 TRUE                     356
3 NA                        20
'
  
  #length     n
  #<chr>  <int>
  #  1 eh      5695
  #2 trgt   25933

#####
##################
#################


library(tidyverse)
#eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% select(STR_ID) %>% count(STR_ID) %>% count(n)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% select(-TRGT_AM) -> eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC

trgt_cov <- read_table("~/Desktop/KU/@research/STR/trgt/oribam_cover/trgt_normal_patho_oribam_coverage.txt") %>% select(-STR_region,-type,-platform) %>% 
  filter(STR_ID %in% simple_STR$STR_ID)
head(trgt_cov)
trgt_cov %>% select(ID) %>% unique() %>% dim()

colnames(trgt_cov) <- c("trgt_meandepth","trgt_meanbaseq","trgt_meanmapq","ID","STR_ID") 
colnames(trgt_cov) <- c("meandepth","meanbaseq","meanmapq","ID","STR_ID") 

eh_cov <- read_table("~/Desktop/KU/@research/STR/eh/eh_simpleSTR_oribam_coverage.txt") %>% select(meandepth:ID) %>% filter(STR_ID %in% simple_STR$STR_ID)
#eh_cov_patho <- read_table("~/Desktop/KU/@research/STR/eh/") #%>% select(meandepth:ID) %>% filter(STR_ID %in% common_STR$STR_DB)

eh_cov %>% select(ID) %>% unique() %>% dim()
head(eh_cov)
colnames(eh_cov) <- c("eh_meandepth","eh_meanbaseq","eh_meanmapq","STR_ID","ID") 
colnames(eh_cov) <- c("meandepth","meanbaseq","meanmapq","STR_ID","ID") 


head(eh_cov)
head(trgt_cov)

eh_cov %>% group_by(STR_ID) %>% summarise(eh_meandepth=mean(meandepth),eh_meanbaseq=mean(meanbaseq),eh_meanmapq=mean(meanmapq)) %>% mutate(platfrom = "eh") -> eh_cov_mean
trgt_cov %>% group_by(STR_ID) %>% summarise(trgt_meandepth=mean(meandepth),trgt_meanbaseq=mean(meanbaseq),trgt_meanmapq=mean(meanmapq)) %>% mutate(platfrom = "trgt") -> trgt_cov_mean

head(eh_cov_mean)
head(trgt_cov_mean)
head(concordance_bySTR_simpleSTR.afterchrXQC)
concordance_bySTR_simpleSTR.afterchrXQC %>% pivot_longer(trgt_meandepth:eh_meanmapq) %>% 
  mutate(platform = ifelse(str_detect(name,"trgt"),"LRS","SRS")) %>%
  group_by(name,platform) %>% mutate(name = str_split_fixed(name,"_",2)[,2]) %>% arrange(name) %>%
  summarise(mean(value))

concordance_bySTR_simpleSTR.afterchrXQC %>% filter(concordance_rate == 0) ->concordance_bySTR_simpleSTR.afterchrXQC_0

concordance_bySTR_simpleSTR.afterchrXQC %>% left_join(trgt_cov_mean %>% select(-platfrom)) %>% left_join(eh_cov_mean %>% select(-platfrom)) %>% #head()
  mutate(concordance_rate_range = case_when(
    concordance_rate == 0 ~ "0",
    concordance_rate > 0 & concordance_rate < 0.1 ~ "(0,0.1)",
    concordance_rate >= 0.1 & concordance_rate < 0.2 ~ "[0.1,0.2)",
    concordance_rate >= 0.2 & concordance_rate < 0.3 ~ "[0.2,0.3)",
    concordance_rate >= 0.3 & concordance_rate < 0.4 ~ "[0.3,0.4)",
    concordance_rate >= 0.4 & concordance_rate < 0.5 ~ "[0.4,0.5)",
    concordance_rate >= 0.5 & concordance_rate < 0.6 ~ "[0.5,0.6)",
    concordance_rate >= 0.6 & concordance_rate < 0.7 ~ "[0.6,0.7)",
    concordance_rate >= 0.7 & concordance_rate < 0.8 ~ "[0.7,0.8)",
    concordance_rate >= 0.8 & concordance_rate < 0.9 ~ "[0.8,0.9)",
    concordance_rate >= 0.9 & concordance_rate < 1 ~ "[0.9,1)",
    concordance_rate == 1 ~ "1")) -> concordance_bySTR_simpleSTR.afterchrXQC
head(concordance_bySTR_simpleSTR.afterchrXQC)
concordance_bySTR_simpleSTR.afterchrXQC %>% write.table("~/Desktop/KU/@research/STR/figure/sup.figure/concordance_bySTR_simpleSTR.bamQC.INFO.afterchrXQC",col.names = T,row.names = F,quote = F,sep = "\t")

concordance_bySTR_simpleSTR.afterchrXQC %>%
  pivot_longer(cols = c("trgt_meanmapq","eh_meanmapq")) %>% #head()
  ggplot(aes(x=value,fill=name)) + 
  geom_density() + 
  facet_wrap(~factor(concordance_rate_range,levels= c("0", "(0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)",
                                                     "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)", "[0.7,0.8)", 
                                                     "[0.8,0.9)", "[0.9,1)","1")),scales = 'free_y',nrow = 3)
head(eh_cov)
head(trgt_cov)
eh_cov %>% filter(STR_ID %in% concordance_bySTR_simpleSTR.afterchrXQC_0$STR_ID) %>% mutate(platform = "eh") %>%
  rbind(trgt_cov %>% filter(STR_ID %in% concordance_bySTR_simpleSTR.afterchrXQC_0$STR_ID) %>% mutate(platform = "trgt")) -> trgt_eh_cov_merge_concordance0
head(trgt_eh_cov_merge_concordance0)

trgt_eh_cov_merge_concordance0 %>% count(platform)
trgt_eh_cov_merge_concordance0 %>% ggplot(aes(x=meanbaseq,fill=platform)) +
  geom_density()




