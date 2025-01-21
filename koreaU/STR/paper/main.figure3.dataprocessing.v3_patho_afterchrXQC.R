## main.figure3 .data processing # patho # after chrX , with PRC1 20241028
library(tidyverse)
library(ggExtra)
library(ggside)

head(final_ref)
dim(final_ref)
final_ref %>% mutate(a= ifelse(str_detect(ID,"chr"),"normal","patho")) %>% count(a)
final_ref %>% filter(!str_detect(ID,"chr")) %>% mutate(a = ifelse(str_detect(MOTIFS,","),"simple","complex")) %>% count(a)

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro
head(final_ref)
final_ref_pro

final_ref %>% filter(str_detect(MOTIFS,",")) %>% 
  separate_rows(MOTIFS, sep = ",") %>%
  rbind(final_ref %>% filter(ID == "RFC1")) %>% #filter(ID == "RFC1")
  group_by(chrom, start, end, ID) %>%
  mutate(priority = row_number()) %>%
  ungroup() %>%
  mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX"))) %>%
  arrange(chrom, start) %>%
  rename(STR_ID = ID) %>%
  mutate(new_ID = paste0(STR_ID,"_",MOTIFS)) -> final_ref_complex
final_ref_complex %>% filter(STR_ID == "RFC1")



sex_info <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/Revio.WGS.sex.info.txt") %>% select(Revio,sex) %>% rename(ID = Revio)
head(sex_info)
sex_info %>% filter(sex == "F") -> sex_info_female
sex_info %>% filter(sex != "F") -> sex_info_male

sex_info %>% count(sex)
head(sex_info_female)
head(patho_complex)
head(patho_complex_male)
patho_complex <- read_table("~/Desktop/KU/@research/STR/figure/complexSTR_prep.EH_TRGT_merged.txt")
patho_complex_male <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/chrX/eh_trgt_complexSTR_prep_chrX.txt") %>% select(ID,STR_ID,new_ID,TRGT_STR,EH_STR,AP,AM)

head(patho_complex)
colnames(patho_complex) <- c("ID","STR_ID","new_ID","MC1","MC2","AM1","AM2","AP1","AP2","EH_STR1","EH_STR2")
head(patho_complex_male)
patho_complex %>% count(MC1 <= MC2)
patho_complex %>% count(EH_STR1 <= EH_STR2)


patho_complex %>% mutate(TRGT_STR1 = ifelse(MC1 < MC2, MC1,MC2),TRGT_STR2 = ifelse(MC1 < MC2, MC2,MC1)) %>%
  mutate(TRGT_AM1 = ifelse(MC1 < MC2, AM1,AM2),TRGT_AM2 = ifelse(MC1 < MC2, AM2,AM1)) %>%
  mutate(TRGT_AP1 = ifelse(MC1 < MC2, AP1,AP2),TRGT_AP2 = ifelse(MC1 < MC2, AP2,AP1)) %>% 
  select(-MC1,-MC2,-AM1,-AM2,-AP1,-AP2) -> patho_complex

patho_complex %>% count(TRGT_STR1 <= TRGT_STR2)

patho_complex %>% select(ID,STR_ID,new_ID,TRGT_STR1,EH_STR1,TRGT_AP1,TRGT_AM1) %>% mutate(allele = 1) -> patho_complex_1
patho_complex %>% select(ID,STR_ID,new_ID,TRGT_STR2,EH_STR2,TRGT_AP2,TRGT_AM2) %>% mutate(allele = 2)-> patho_complex_2
colnames(patho_complex_1) <- c("ID","STR_ID","new_ID","TRGT_STR","EH_STR","AP","AM","allele")
colnames(patho_complex_2) <- c("ID","STR_ID","new_ID","TRGT_STR","EH_STR","AP","AM","allele")
patho_complex_1 %>% rbind(patho_complex_2) %>% filter(ID %in% sex_info_male$ID & STR_ID %in% patho_complex_male$STR_ID)

patho_complex_1 %>% rbind(patho_complex_2) %>% filter(!(ID %in% sex_info_male$ID & STR_ID %in% patho_complex_male$STR_ID)) %>% #count(ID) %>% count(n)
  rbind(patho_complex_male %>% mutate(allele = 1)) -> patho_complex_processing

patho_complex_processing %>% #filter(STR_ID == "HTT") %>%
  mutate(TRGT_STR = ifelse(new_ID == "HTT_CAG",TRGT_STR-2,EH_STR)) -> patho_complex_processing

head(patho_complex_processing)

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.txt")

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1<- eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1")

add_RFC_TRGT <- read_table("~/Desktop/KU/@research/STR/trgt/Quality_check_RFC1/RFC_merge_pro.txt")
head(add_RFC_TRGT)
add_RFC_TRGT %>% count(TRGT_STR1 <=TRGT_STR2)

add_RFC_TRGT %>% filter(ID %in% patho_complex_processing$ID) %>% select(ID,TRGT_STR1,TRGT_length1,TRGT_AP1) %>% mutate(allele = 1) -> add_RFC_TRGT1
add_RFC_TRGT %>% filter(ID %in% patho_complex_processing$ID) %>% select(ID,TRGT_STR2,TRGT_length2,TRGT_AP2) %>% mutate(allele = 2) -> add_RFC_TRGT2

colnames(add_RFC_TRGT1) <- c("ID","TRGT_STR","TRGT_length","TRGT_AP","allele")
colnames(add_RFC_TRGT2) <- c("ID","TRGT_STR","TRGT_length","TRGT_AP","allele")

head(patho_complex_1)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1)
add_RFC_TRGT1 %>% rbind(add_RFC_TRGT2) %>% mutate(STR_ID ="RFC1",new_ID = "RFC1_AARRG") %>% mutate(AM = "-1") %>%
  left_join(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1 %>% select(ID,STR_ID,EH_STR,allele)) -> eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1_reQC
#head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1_reQC)
#write.table(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1_reQC,"")


head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID == "RFC1")

head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1)
head(patho_complex_processing)



eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(str_detect(STR_ID,"chr")) %>% 
  filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(STR_ID) %>% 
  summarise(concordance_rate = mean(check)) -> concordance_bySTR_simpleSTR_normalSTR

head(concordance_bySTR_simpleSTR_normalSTR)

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(str_detect(STR_ID,"chr")) %>%
  filter(STR_ID != "RFC1") %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID) %>% 
  summarise(concordance_rate = mean(check)) -> concordance_byID_simpleSTR_normalSTR

head(concordance_byID_simpleSTR_normalSTR)

head(concordance_byID_simpleSTR_normalSTR)
head(concordance_bySTR_simpleSTR_normalSTR)
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)
eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>% filter(STR_ID != "RFC1") %>%
  filter(!str_detect(STR_ID,"chr")) -> eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_onlysimplepath

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_onlysimplepath %>% write.table("~/Desktop/KU/@research/STR/figure/simplepathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_onlysimplepath)

head(patho_simple)

#patho_simple %>% filter(STR_ID == "RFC1")
patho_simple <- read_table("~/Desktop/KU/@research/STR/figure/simplepathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")
patho_complex_processing$STR_type <- "complex"
patho_simple$STR_type <- "simple"

head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1)
head(patho_complex_processing)


head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1_reQC)
head(patho_complex_processing)

#patho_complex_processing %>% #head()
#  rbind(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1 %>% mutate(new_ID = STR_ID) %>% rename(AM = TRGT_AM,AP = TRGT_AP) %>% mutate(STR_type = "complex")) %>%
#  write.table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")


# RFC1
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1_reQC)
head(patho_complex_processing)
patho_complex_processing %>% #head()
  rbind(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC_RFC1_reQC %>%  mutate(STR_type = "complex") %>%
          select(-TRGT_length) %>% rename(AP = TRGT_AP)) %>% 
  write.table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")


head(patho_complex_processing)
patho_complex_processing %>% filter(STR_ID == "RFC1")
head(patho_simple)
patho_complex_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")
#head(patho_complex_processing)
patho_complex_processing %>% filter(STR_ID == "RFC1")
patho_complex_processing %>% filter(TRGT_STR != EH_STR) -> a
patho_complex_processing %>% filter(STR_ID == "CNBP") %>% filter(ID == "NIH23F1777753")#filter(ID == "NIH20N2552908")
head(patho_complex_processing)
patho_complex_processing %>%  select(STR_ID) %>% unique()

patho_simple %>% mutate(new_ID = STR_ID) %>% rename(AM = TRGT_AM,AP = TRGT_AP)  %>% rbind(patho_complex_processing) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(new_ID) %>% 
  summarise(concordance_rate = mean(check)) -> a

#NIH23F1110753 ATXN8OS
patho_simple %>% filter(STR_ID == "RFC1")
head(patho_complex_processing)
patho_complex_processing
patho_complex_processing %>% filter(TRGT_STR != EH_STR) -> a

patho_simple %>% mutate(new_ID = STR_ID) %>% rename(AM = TRGT_AM,AP = TRGT_AP)  %>% rbind(patho_complex_processing) %>% #filter(ID == "NIH23F1110753" & STR_ID == "ATXN8OS")
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID,STR_ID,allele) %>% #filter(STR_ID != "RFC1") %>%
  summarise(check = mean(check)) %>% #filter(check %in% c(0.5))
  mutate(check = ifelse(check == 1,1,0)) %>% group_by(ID) %>% #head()
  summarise(concordance_rate = mean(check)) -> concordance_byID_pathoSTR

patho_simple %>% mutate(new_ID = STR_ID) %>% rename(AM = TRGT_AM,AP = TRGT_AP)  %>% rbind(patho_complex_processing) %>% #filter(ID == "NIH23F1110753" & STR_ID == "ATXN8OS")
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID,STR_ID,allele) %>% #filter(STR_ID != "RFC1") %>%
  summarise(check = mean(check)) %>% #filter(check %in% c(0.5))
  mutate(check = ifelse(check == 1,1,0)) %>% group_by(STR_ID) %>% #head()
  summarise(concordance_rate = mean(check)) -> concordance_bySTR_pathoSTR


#concordance_byID_complexSTR_pathoSTR

patho_simple %>% mutate(new_ID = STR_ID) %>% rename(AM = TRGT_AM,AP = TRGT_AP)  %>% #rbind(patho_complex_processing) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID) %>% filter(STR_ID != "RFC1") %>%
  summarise(concordance_rate = mean(check)) -> concordance_byID_pathoSTR_simple
#concordance_byID_pathoSTR_simple

patho_complex_processing %>% mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID,STR_ID,allele) %>% #filter(STR_ID != "RFC1") %>%
  summarise(check = mean(check)) %>% #filter(check %in% c(0.5))
  mutate(check = ifelse(check == 1,1,0)) %>% group_by(ID) %>% #head()
  summarise(concordance_rate = mean(check)) -> concordance_byID_pathoSTR_complex


patho_simple %>% mutate(new_ID = STR_ID) %>% rename(AM = TRGT_AM,AP = TRGT_AP)  %>% #rbind(patho_complex_processing) %>% #head()
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(STR_ID) %>% filter(STR_ID != "RFC1") %>%
  summarise(concordance_rate = mean(check)) -> concordance_bySTR_pathoSTR_simple

head(patho_complex_processing)
patho_complex_processing %>% filter(STR_ID == "RFC1") %>% filter(TRGT_STR !=EH_STR)
patho_complex_processing %>% filter(TRGT_STR !=EH_STR) -> a

patho_complex_processing %>% mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% group_by(ID,STR_ID,allele) %>% #filter(STR_ID != "RFC1") %>%
  summarise(check = mean(check)) %>% #filter(check %in% c(0.5))
  mutate(check = ifelse(check == 1,1,0)) %>% group_by(STR_ID) %>% #head()
  summarise(concordance_rate = mean(check)) -> concordance_bySTR_pathoSTR_complex



head(concordance_byID_pathoSTR_simple)
head(concordance_byID_pathoSTR_complex)
head(concordance_byID_simpleSTR_normalSTR)
head(concordance_byID_pathoSTR)

concordance_byID_pathoSTR_simple$g = "Pathogenic"
concordance_byID_pathoSTR_complex$g = "Pathogenic"
concordance_byID_pathoSTR$g = "Pathogenic"
concordance_byID_simpleSTR_normalSTR$g = "Normal"

concordance_byID_pathoSTR_complex$str_type = "(complexSTR)"
concordance_byID_pathoSTR_simple$str_type = "(simpleSTR)"
concordance_byID_pathoSTR$str_type = "(Overall)"
concordance_byID_simpleSTR_normalSTR$str_type = "(simpleSTR)"

head(concordance_byID_pathoSTR_complex)
concordance_byID_pathoSTR_complex %>% 
  rbind(concordance_byID_pathoSTR_simple) %>%
  rbind(concordance_byID_pathoSTR) %>%
  rbind(concordance_byID_simpleSTR_normalSTR) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure3/concordance.byID.nonpatho.patho.simpleSTR.complexSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  mutate(new_g = paste0(g,"\n",str_type)) %>%
  ggplot(aes(x=factor(new_g,levels=c("Normal\n(simpleSTR)","Pathogenic\n(Overall)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")),y=concordance_rate,fill=factor(new_g,levels=c("non-Pathogenic\n(simpleSTR)","Pathogenic\n(Overall)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")))) + 
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black",alpha=0.8) + 
  theme_bw() +
#  theme_step1() + 
  theme(legend.position = 'none',
        axis.title.x=element_blank())


#####
head(concordance_bySTR_simpleSTR_normalSTR)
head(concordance_bySTR_pathoSTR)
head(concordance_bySTR_pathoSTR_simple)
head(concordance_bySTR_pathoSTR_complex)

concordance_bySTR_pathoSTR_simple$g = "Pathogenic"
concordance_bySTR_pathoSTR_complex$g = "Pathogenic"
concordance_bySTR_pathoSTR$g = "Pathogenic"
concordance_bySTR_simpleSTR_normalSTR$g = "Normal"

concordance_bySTR_pathoSTR_complex$str_type = "(complexSTR)"
concordance_bySTR_pathoSTR_simple$str_type = "(simpleSTR)"
concordance_bySTR_pathoSTR$str_type = "(Overall)"
concordance_bySTR_simpleSTR_normalSTR$str_type = "(simpleSTR)"


concordance_bySTR_pathoSTR_complex %>% 
  rbind(concordance_bySTR_pathoSTR_simple) %>%
  rbind(concordance_bySTR_pathoSTR) %>%
  rbind(concordance_bySTR_simpleSTR_normalSTR) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure3/concordance.bySTR.nonpatho.patho.simpleSTR.complexSTR.afterchrXQC.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  mutate(new_g = paste0(g,"\n",str_type)) %>%
  ggplot(aes(x=factor(new_g,levels=c("Normal\n(simpleSTR)","Pathogenic\n(Overall)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")),y=concordance_rate,fill=factor(new_g,levels=c("Normal\n(simpleSTR)","Pathogenic\n(Overall)","Pathogenic\n(simpleSTR)","Pathogenic\n(complexSTR)")))) + 
  geom_boxplot() + 
#  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, color = "black",alpha=0.8) + 
  theme_bw() +
  #  theme_step1() + 
  theme(legend.position = 'none',
        axis.title.x=element_blank())


###

head(final_ref_complex)
final_ref_complex
patho_complex_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")
head(patho_complex_processing)
patho_complex_processing %>% filter(STR_ID == "RFC1")
final_ref_complex
patho_complex_processing %>% #head()
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>% #head()
  group_by(ID,STR_ID,allele) %>%
  summarise(check = mean(check)) %>% mutate(check = ifelse(check == 1,1,0)) %>% ungroup() %>%
  count(STR_ID,check) %>% group_by(STR_ID) %>%
  mutate(prop = prop.table(n)) %>% #group_by(STR_ID,check) %>%
  mutate(ref = if_else(check == 0, n, NA_real_)) %>%
  fill(ref, .direction = "downup") %>%
  mutate(ref = ifelse(is.na(ref),0,ref))  %>%
#  mutate(STR_ID = fct_reorder(STR_ID, ref)) %>%
  ggplot(aes(x = fct_reorder(STR_ID,ref), y = n, fill = factor(check))) + 
  geom_bar(stat = 'identity', position = "fill") + 
  theme_bw() +
  xlab("STR_ID (sorted by ref)") +
  ylab("Proportion")
    
head(final_ref_complex)
head(patho_complex_processing)
head(patho_simple)
patho_simple %>% left_join(final_ref_pro) %>%
  mutate(TRGT_STR_length = TRGT_STR*str_length(MOTIFS),EH_STR_length = EH_STR*str_length(MOTIFS)) %>% #head()
  mutate(check= ifelse(TRGT_STR_length == EH_STR_length, 1, 0)) %>% 
  select(ID,STR_ID,TRGT_STR_length,EH_STR_length,TRGT_AP,allele,STR_type,check)
    
  
  
head(patho_complex_processing)
patho_complex_processing %>% filter(TRGT_STR < EH_STR)
patho_complex_processing %>% filter(STR_ID == "RFC1")
head(patho_simple)
patho_simple %>% filter(TRGT_STR - EH_STR < -10)
patho_complex_processing %>% #filter(STR_ID == "RFC1")
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>% #filter(STR_ID == 'RFC1') #%>%
  left_join(final_ref_complex %>% select(STR_ID,new_ID,MOTIFS)) %>% #filter(STR_ID == "RFC1")
  mutate(TRGT_STR_length = TRGT_STR*str_length(MOTIFS),EH_STR_length = EH_STR*str_length(MOTIFS)) %>%
  group_by(ID,STR_ID,allele) %>% #filter(AP < 0)
  summarise(TRGT_STR_length = sum(TRGT_STR_length),EH_STR_length = sum(EH_STR_length),TRGT_AP = mean(AP)) %>% 
  mutate(check= ifelse(TRGT_STR_length == EH_STR_length, 1, 0),STR_type="complex") %>% #head()
  rbind(patho_simple %>% left_join(final_ref_pro) %>%
          mutate(TRGT_STR_length = TRGT_STR*str_length(MOTIFS),EH_STR_length = EH_STR*str_length(MOTIFS)) %>% #head()
          mutate(check= ifelse(TRGT_STR_length == EH_STR_length, 1, 0)) %>% 
          select(ID,STR_ID,TRGT_STR_length,EH_STR_length,TRGT_AP,allele,STR_type,check)) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure3/pathogenic.all.allele.strlength.compare.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  ggplot(aes(x=TRGT_STR_length,y=EH_STR_length,shape=STR_type,color = factor(check))) +
  geom_point(alpha=0.7) +
  geom_xsidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  
  geom_ysidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) + 
  theme_bw() +
  coord_equal(ratio = 1)
   

patho_complex_processing %>% 
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>% #filter(STR_ID != 'RFC1') %>%
  group_by(STR_ID, new_ID) %>%
  summarise(concordance_rate = mean(check)) %>% 
  left_join(final_ref_complex %>%   mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX")))) %>% # head()
  group_by(STR_ID) %>% 
  mutate(MOTIFS = factor(MOTIFS, levels = unique(final_ref_complex$MOTIFS[order(final_ref_complex$priority)]))) %>% 
  filter(STR_ID != "PRNP") %>%
  ggplot(aes(x = MOTIFS, y = concordance_rate, fill = STR_ID)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~factor(STR_ID,levels = c("ATXN7","HTT","FXN","FRA10AC1","ATXN8OS","NOP56",#"PRNP",
                                       "TMEM185A","CNBP")),scales = "free_x",ncol = 4) + 
  theme_minimal() 


head(patho_complex_processing_each)
head(patho_complex_processing_each)
patho_complex_processing %>% 
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>%  #head()
  group_by(ID,STR_ID,allele) %>% 
  summarise(concordance_rate = mean(check)) %>% 
  mutate(check= ifelse(concordance_rate == 1, 1, 0)) %>%  
  group_by(STR_ID) %>%  
  summarise(concordance_rate = mean(check)) -> patho_complex_processing_overall

head(patho_complex_processing_overall)
patho_complex_processing_overall
patho_complex_processing %>% 
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>%
  group_by(STR_ID, new_ID) %>%
  summarise(concordance_rate = mean(check)) %>% 
  left_join(final_ref_complex %>%   mutate(chrom = factor(chrom, levels = c(paste0("chr", 1:22), "chrX")))) %>% # head()
  group_by(STR_ID) %>% 
  mutate(MOTIFS = factor(MOTIFS, levels = unique(final_ref_complex$MOTIFS[order(final_ref_complex$priority)]))) %>% 
  filter(STR_ID != "PRNP") %>%
  left_join(patho_complex_processing_overall %>% rename(mean_line = concordance_rate)) %>%  
  ggplot(aes(x = MOTIFS, y = concordance_rate, fill = STR_ID)) + 
  geom_bar(stat = "identity") + 
  geom_hline(aes(yintercept = mean_line), linetype = "dashed", color = "red") +  
  facet_wrap(~factor(STR_ID,levels = c("ATXN7","HTT","FXN","FRA10AC1","ATXN8OS","NOP56",#"PRNP",
                                       "TMEM185A","CNBP")),scales = "free_x",ncol = 4) + 
  theme_minimal() 




## suple
patho_complex_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")

patho_complex_processing %>% 
  mutate(check= ifelse(TRGT_STR != EH_STR, 1, 0)) %>% 
  #count(ID,STR_ID)
  group_by(ID,STR_ID,allele) %>%
  summarise(check_rate = sum(check)) %>% #filter(STR_ID == "CNBP")
  group_by(STR_ID) %>% filter(STR_ID != "PRNP") %>%
  count(check_rate) %>% #head()
  ggplot(aes(x=factor(STR_ID,levels = c("ATXN7","HTT","FXN","FRA10AC1","ATXN8OS","NOP56",#"PRNP",
                                        "TMEM185A","CNBP")),y=n,fill=rev(factor(check_rate)))) + 
  geom_bar(stat = 'identity',position = "fill")

  
patho_complex_processing %>% #head()
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>% 
  #count(ID,STR_ID)
  group_by(ID,STR_ID,allele) %>%
  summarise(check_rate = sum(check)) %>% #filter(STR_ID == "CNBP")
  group_by(STR_ID) %>% #filter(STR_ID != "PRNP") %>% 
  count(check_rate) %>% mutate(prop = prop.table(n)*100) %>% mutate(n_prop = paste0(n," (",round(prop,1),"%)")) %>%
  arrange(check_rate) %>% select(-n,-prop) %>%
  pivot_wider(names_from = check_rate,values_from = n_prop) %>% 
  writexl::write_xlsx("~/Desktop/KU/@research/STR/table/complex.STR.each.motif.match.count.xlsx")
  

#### patho mapQ

patho_simple <- read_table("~/Desktop/KU/@research/STR/figure/simplepathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt") #%>% filter(STR_ID == "RFC1")
head(patho_simple)
patho_simple %>% filter(STR_ID == "RFC1")
patho_complex_processing <- read_table("~/Desktop/KU/@research/STR/figure/figure3/complex_pathogenic_STR_prep.EH_TRGT_merged_afterchrXQC.txt")
head(patho_simple)

patho_complex_processing %>% filter(STR_ID == "RFC1")
head(patho_complex_processing)
patho_complex_processing %>% select(STR_ID) %>% unique()

trgt_eh_patho_oribam <- read_table("~/Desktop/KU/@research/STR/figure/qc/trgt_eh_patho_oribam.txt")
head(trgt_eh_patho_oribam)

trgt_eh_patho_oribam %>% pivot_longer(meandepth:meanmapq) %>% #head()
  mutate(platform = ifelse(platform == "LRS","TRGT","EH")) %>%
  pivot_wider(names_from = c(platform, name), values_from = value) 

patho_complex_processing %>% left_join(final_ref_complex %>% select(STR_ID,new_ID,MOTIFS)) %>% mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(TRGT_STR_length = RU.length*TRGT_STR,EH_STR_length=RU.length*EH_STR) %>%
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>% group_by(ID,STR_ID,allele) %>% 
  summarise(TRGT_AP = mean(AP),check = mean(check),TRGT_STR_length_sum= sum(TRGT_STR_length),EH_STR_length_sum= sum(EH_STR_length)) %>%
  mutate(check=ifelse(check == 1,1,0)) %>% mutate(type = "complex") -> patho_complex_processing


head(patho_complex_processing)
head(patho_simple)

patho_simple %>% left_join(final_ref_pro) %>% mutate(TRGT_STR_length = RU.length*TRGT_STR,EH_STR_length=RU.length*EH_STR) %>%
  mutate(check= ifelse(TRGT_STR == EH_STR, 1, 0)) %>% group_by(ID,STR_ID,allele) %>% 
  summarise(TRGT_AP = mean(TRGT_AP),check = mean(check),TRGT_STR_length_sum= sum(TRGT_STR_length),EH_STR_length_sum= sum(EH_STR_length)) %>% #head()
  mutate(type = "simple") -> patho_simple

head(patho_simple)
patho_simple %>% rbind(patho_complex_processing) %>% 
  rename(TRGT_STR_length = TRGT_STR_length_sum,EH_STR_length = EH_STR_length_sum) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure4/pathoSTR.STRlength.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  left_join(trgt_eh_patho_oribam %>% pivot_longer(meandepth:meanmapq) %>% #head()
              mutate(platform = ifelse(platform == "LRS","trgt","eh")) %>%
              pivot_wider(names_from = c(platform, name), values_from = value)) %>% write.table("~/Desktop/KU/@research/STR/figure/figure4/pathoSTR.STRlength.with.bamstats.txt",col.names = T,row.names = F,quote = F,sep = "\t")

  
patho_simple %>% rbind(patho_complex_processing) %>% 
  left_join(trgt_eh_patho_oribam %>% select(meanbaseq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanbaseq)) %>%
  ggplot(aes(x=LRS,y=SRS,color=type)) +
  geom_point(alpha=0.8) + 
  facet_grid(~check)
head(trgt_eh_patho_oribam)

patho_simple %>% rbind(patho_complex_processing) %>% #head()
  left_join(trgt_eh_patho_oribam %>% select(meanmapq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanmapq)) %>%
  ggplot(aes(x=LRS,y=SRS,color=type)) +
  geom_point()

patho_simple %>% rbind(patho_complex_processing) %>% #head()
  left_join(trgt_eh_patho_oribam %>% select(meanmapq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanmapq)) %>%
  ggplot(aes(x=LRS,y=SRS,color=type)) +
  geom_point() +
  facet_grid(~check)

patho_simple %>% rbind(patho_complex_processing) %>% #head()
  left_join(trgt_eh_patho_oribam %>% select(meanmapq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanmapq)) %>%
  group_by(type,check,STR_ID) %>% 
  summarise(LRS = mean(LRS), SRS = mean(SRS)) %>%
  ggplot(aes(x=LRS,y=SRS)) +
  #geom_point() + 
  geom_point(aes(col = factor(check),shape=type)) +  # Scatter plot with colors by 'chec
  #geom_point(aes(shape = type)) +  # Scatter plot with colors by 'check'
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Red dashed y = x line
  #xlim(c(25,40)) +  # Set x-axis limits
  #ylim(c(25,40)) + 
  geom_xsidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on x-axis
  geom_ysidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on y-axis
  theme_bw() + 
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.panel.scale = .1
  )



patho_simple %>% rbind(patho_complex_processing) %>% #head()
  left_join(trgt_eh_patho_oribam %>% select(meanmapq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanmapq)) %>%
  ggplot(aes(x=LRS,y=SRS)) +
  #geom_point() + 
  geom_point(aes(col = factor(check),shape=type)) +  # Scatter plot with colors by 'chec
  #geom_point(aes(shape = type)) +  # Scatter plot with colors by 'check'
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Red dashed y = x line
  #xlim(c(25,40)) +  # Set x-axis limits
  #ylim(c(25,40)) + 
  geom_xsidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on x-axis
  geom_ysidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on y-axis
  theme_bw() + 
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.panel.scale = .1
  )





patho_complex_processing %>% #head()
  left_join(trgt_eh_patho_oribam %>% select(meanbaseq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanbaseq)) %>%
  ggplot(aes(x=LRS,y=SRS,color=check)) +
  geom_point()




library(palmerpenguins)
library(ggside)

patho_simple %>% rbind(patho_complex_processing) %>% 
  left_join(trgt_eh_patho_oribam %>% select(meanbaseq,ID,STR_ID,platform) %>% 
              pivot_wider(names_from = platform,values_from = meanbaseq)) %>% 
  ggplot(aes(x=LRS,y=SRS)) +
  geom_point(aes(col=factor(check))) + 
  #geom_point() + 
  geom_xsidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) + 
  theme_bw()
## good
patho_simple %>% rbind(patho_complex_processing) %>%   
  left_join(trgt_eh_patho_oribam %>%
              select(meanbaseq, ID, STR_ID, platform) %>%
              pivot_wider(names_from = platform, values_from = meanbaseq)) %>% #write.table("~/Desktop/KU/@research/STR/figure/figure3/meanbaseQ.patho_simple_complex.all.allele.txt",col.names = T,row.names = F,quote = F,sep = "\t")
  ggplot(aes(x = LRS, y = SRS)) +
  geom_point(aes(col = factor(check))) +  # Scatter plot with colors by 'chec
  #geom_point(aes(shape = type)) +  # Scatter plot with colors by 'check'
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Red dashed y = x line
  xlim(c(25,40)) +  # Set x-axis limits
  ylim(c(25,40)) + 
  geom_xsidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on x-axis
  geom_ysidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on y-axis
  theme_bw() + 
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.panel.scale = .1
  )

patho_simple %>% #rbind(patho_complex_processing) %>%   
  left_join(trgt_eh_patho_oribam %>%
              select(meanbaseq, ID, STR_ID, platform) %>%
              pivot_wider(names_from = platform, values_from = meanbaseq)) %>% #head()
  group_by(type,check,STR_ID) %>% 
  summarise(LRS = mean(LRS), SRS = mean(SRS)) %>%
  ggplot(aes(x = LRS, y = SRS)) +
  geom_point(aes(col = factor(check))) +  # Scatter plot with colors by 'chec
  #geom_point(aes(shape = type)) +  # Scatter plot with colors by 'check'
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Red dashed y = x line
  xlim(c(25,40)) +  # Set x-axis limits
  ylim(c(25,40)) + 
  geom_xsidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on x-axis
  geom_ysidedensity(aes(fill = factor(check)), alpha = 0.5, show.legend = FALSE) +  # Density on y-axis
  theme_bw() + 
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.panel.scale = .1
  )





patho_complex_processing %>% #head()
  left_join(trgt_eh_patho_oribam %>% select(meanmapq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanmapq)) %>%
  ggplot(aes(x=LRS,y=SRS,color=check)) +
  geom_point()
  
patho_complex_processing %>% #head()
  left_join(trgt_eh_patho_oribam %>% select(meanmapq,ID,STR_ID,platform) %>% pivot_wider(names_from = platform,values_from = meanmapq)) %>% #head()
    ggplot(aes(x=TRGT_STR_length_sum,y=EH_STR_length_sum,color=TRGT_AP))  +
    geom_point()

