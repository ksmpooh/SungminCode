## table
library(tidyverse)
library(writexl)
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/table")
## concordacne rate

final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

final_ref %>% select(MOTIFS,ID) %>% rename(STR_ID = ID) %>% 
  mutate(RU.length = str_length(MOTIFS)) %>%
  mutate(GC = str_length(gsub("[^CG]", "",MOTIFS))/str_length(MOTIFS)) -> final_ref_pro

ref <- read.table("~/Desktop/KU/@research/STR/db/annovar/output_bed_ref")
ref %>% unique()-> ref

colnames(ref) <- c("chrom","start","end","STR_ID")


#simple_concordance <- read.table("~/Desktop/KU/@research/STR/figure/figure2/f2.concordance_bySTR_simpleSTR.txt",header = T)
simple_concordance <- read.table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.concordance_bySTR_simpleSTR.afterchrXQC.txt",header = T)
#head(simple_concordance)
apscore <- read_table("~/Desktop/KU/@research/STR/trgt/trgt_simpleSTR_intersect.AP.score.byAllele.txt")
head(apscore)
apscore %>% filter(STR_ID != "RFC1") %>%
  group_by(STR_ID) %>% summarise(APmean = mean(AP)) -> apscore_mean
head(simple_concordance)
mean(simple_concordance$concordance_rate) #0.9239236
#simple_concordance %>% mutate(g = ifelse(str_detect(STR_ID,"chr"),"normal","patho")) %>% group_by(g) %>% #head()

simple_concordance %>% left_join(apscore_mean) %>% mutate(pathogenic= ifelse(str_detect(STR_ID,"chr"),"normal","pathogenic")) %>%
  mutate(AP0.5 = ifelse(APmean >=0.5,">=0.5","<0.5")) %>% 
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
    concordance_rate == 1 ~ "1")) %>% count(concordance_rate_range,AP0.5,pathogenic) %>% 
  pivot_wider(names_from = pathogenic,values_from = n) %>%  mutate(across(everything(), ~replace_na(., 0))) %>% 
  pivot_wider(names_from = AP0.5,values_from = c(normal,pathogenic)) %>%
  write_xlsx("simpleSTR.count.byconcordanceRange.pathogenic.APscore.xlsx")


head(ref)
## anno table
anno <- read_table("~/Desktop/KU/@research/STR/db/annovar/Raw.anno.merge.processing.onlyneed.txt") %>% filter(STR_ID != "RFC1")
anno %>% mutate(arm = ifelse(str_detect(cytoband,"p"),"p","q")) %>% 
  mutate(type = ifelse(type == "upstream","promoter",type)) %>% left_join(simple_concordance) %>% na.omit() %>% left_join(ref) %>%
  rename(chromosome = chrom)-> anno.merge
  #mutate(type = ifelse(str_detect(type,"RNA"),"ncRNA",type))
head(anno.merge)
dim(anno.merge)

anno.merge %>% filter(STR_ID != "RFC1") %>%
  group_by(type) %>% summarise(concordance = mean(concordance_rate))
anno.merge %>% 
  count(type) %>% left_join(anno.merge %>% group_by(type) %>% summarise(concordance = mean(concordance_rate))) %>% rename(Annotation = type) -> anno.concordance
anno.merge %>% count(chromosome) %>% left_join(anno.merge %>% group_by(chromosome) %>% summarise(concordance = mean(concordance_rate))) %>% 
  mutate(chromosome = factor(chromosome, levels = c(paste0("chr", 1:22), "chrX"))) %>%
  arrange(chromosome) -> chrom.concordance
anno.merge %>% count(chromosome,arm) %>% left_join(anno.merge %>% group_by(chromosome,arm) %>% summarise(concordance = mean(concordance_rate))) %>%
  mutate(chromosome = factor(chromosome, levels = c(paste0("chr", 1:22), "chrX"))) %>%
  arrange(chromosome)-> chrom.arm.concordance

anno.merge %>% count(arm) %>% left_join(anno.merge %>% group_by(arm) %>% summarise(concordance = mean(concordance_rate))) -> arm.concordance


anno.merge %>% count(main_cyto_type) %>% left_join(anno.merge %>% group_by(main_cyto_type) %>% summarise(concordance = mean(concordance_rate))) %>%
  mutate(main_cyto_type = case_when(
    main_cyto_type == "acen" ~ "Centromere",
    main_cyto_type == "gneg" ~ "G-negative",
    main_cyto_type == "gpos" ~ "G-positive",
    main_cyto_type == "gvar" ~ "G-variable",
    main_cyto_type == "stalk" ~ "Stalk",
    TRUE ~ main_cyto_type  # 해당하지 않는 경우 기존 값 유지
  )) %>% rename(Cytoband = main_cyto_type)  -> cyto.concordance

anno.merge %>% filter(main_cyto_type == "gpos") %>% #count(cyto_type)
  count(cyto_type) %>% left_join(anno.merge %>% filter(main_cyto_type == "gpos") %>% group_by(cyto_type) %>% summarise(concordance = mean(concordance_rate))) %>% 
  mutate(cyto_type = case_when(
    cyto_type == "gpos100" ~ "G-positive 100%",
    cyto_type == "gpos25" ~ "G-positive 25%",
    cyto_type == "gpos50" ~ "G-positive 50%",
    cyto_type == "gpos75" ~ "G-positive 75%",
    TRUE ~ cyto_type  # 해당하지 않는 경우 기존 값 유지
  )) %>%
    mutate(cyto_type = factor(cyto_type, levels = c("G-positive 100%","G-positive 75%","G-positive 50%","G-positive 25%"))) %>%
    arrange(cyto_type) %>% 
    rename(Cytoband_gpos = cyto_type)  -> cyto.gpos.concordance

anno.concordance
chrom.concordance
chrom.arm.concordance
cyto.concordance
cyto.gpos.concordance
#acen: Acentric region (Centromere) 
#gneg: G-negative (Giemsa negative)
#gpos: G-positive (Giemsa positive)
#gvar: G-variable
#stalk: Stalk region 

#acen: Centromere)
#gneg: G-negative 
#gpos: G-positive
#gvar: G-variable
#stalk: Stalk
tail(chrom.concordance)
sheets <- list(
  anno.concordance,
  chrom.concordance,
  arm.concordance,
  chrom.arm.concordance,
  cyto.concordance,
  cyto.gpos.concordance
)

write_xlsx(sheets, path = "simpleSTR.annotation.concordance.xlsx")


anno.concordance$type = "Annotation"
chrom.concordance$type = "Chromosome"
arm.concordance$type = "Arm"
chrom.arm.concordance$type = "Chromosome.Arm"
cyto.concordance$type = "Cytoband"
cyto.gpos.concordance$type = "Cytoband.gpos"
head(chrom.arm.concordance)
colnames(anno.concordance)[1] <- "name"
colnames(chrom.concordance)[1] <- "name"
colnames(arm.concordance)[1] <- "name"
colnames(chrom.arm.concordance)[1] <- "name"
colnames(cyto.concordance)[1] <- "name"
colnames(cyto.gpos.concordance)[1] <- "name"
cyto.gpos.concordance %>% mutate(name = str_replace_all(name,"G-positive ",""))
anno.concordance %>% rbind(chrom.concordance) %>%
  rbind(arm.concordance) %>%
  rbind(cyto.concordance) %>%
  rbind(cyto.gpos.concordance %>% mutate(name = str_replace_all(name,"G-positive ",""))) %>% mutate(arm = "NA") %>%
  rbind(chrom.arm.concordance) %>% write.table("~/Desktop/KU/@research/STR/figure/figure2_withchrX/f2.annotation.info.concordance.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#######   allele count

library(openxlsx)
wb <- createWorkbook()

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/df_merge/eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC.withoutRFC1.txt")
head(eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC)

eh_trgt_merge_simple_pass_intersect_forConcordance.afterchrXQC %>%
  mutate(check = ifelse(TRGT_STR == EH_STR,1,0)) %>% left_join(final_ref_pro) %>% #head()
  mutate(STR_length = TRGT_STR * RU.length) -> basic.count.STRallele
  


head(simpleSTR_count_bylength)
basic.count.STRallele %>% mutate(STR_length = case_when(
  STR_length < 50 ~ "[0~50)",
  STR_length < 100 ~ "[50~100)",
  STR_length < 150 ~ "[100~150)",
  TRUE ~ '[150~Inf)')) %>% count(check,STR_length) %>% 
  pivot_wider(names_from = STR_length,values_from = n) %>% 
  mutate(overall = rowSums(select(., `[0~50)`, `[50~100)`, `[100~150)`, `[150~Inf)`)))

  

head(simpleSTR_count_bylength)
basic.count.STRallele %>% mutate(STR_length = case_when(
  STR_length < 50 ~ "[0~50)",
  STR_length < 100 ~ "[50~100)",
  STR_length < 150 ~ "[100~150)",
  TRUE ~ '[150~Inf)')) %>% count(check,STR_length) %>% 
  pivot_wider(names_from = STR_length,values_from = n) %>% 
  rename(match = check) %>%
  select(match,"[0~50)","[50~100)","[100~150)",'[150~Inf)') -> a
head(a)

wb <- createWorkbook()

addWorksheet(wb, "simpleSTR_basic_count")
writeData(wb, "simpleSTR_basic_count", a)

basic.count.STRallele %>% 
  mutate(GC = case_when(
    GC < 0.25 ~ "[0~0.25)",
    GC < 0.5 ~ "[0.25~0.5)",
    GC < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>% count(check,GC) %>%
  pivot_wider(names_from = GC,values_from = n) %>%
  rename(match = check) %>%
  select(match,"[0~0.25)","[0.25~0.5)","[0.5~0.75)",'[0.75~1]') -> b

addWorksheet(wb, "simpleSTR_GC")
writeData(wb, "simpleSTR_GC", b)


basic.count.STRallele %>% mutate(STR_length = case_when(
  STR_length < 50 ~ "[0~50)",
  STR_length < 100 ~ "[50~100)",
  STR_length < 150 ~ "[100~150)",
  TRUE ~ '[150~Inf)')) %>% 
  mutate(GC = case_when(
    GC < 0.25 ~ "[0~0.25)",
    GC < 0.5 ~ "[0.25~0.5)",
    GC < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>%  count(check,STR_length,GC) -> c
  
  
head(c)

c %>% 
  pivot_wider(names_from = STR_length,values_from = n) %>%
  rename(match = check) %>% mutate(GC = factor(GC,levels = c("[0~0.25)", "[0.25~0.5)", "[0.5~0.75)", "[0.75~1]"))) %>%
  arrange(match, GC) %>%
  select(match,GC,"[0~50)","[50~100)","[100~150)",'[150~Inf)') -> c

addWorksheet(wb, "simpleSTR_STRlengthGC")
writeData(wb, "simpleSTR_STRlengthGC", c)



basic.count.STRallele %>%
  mutate(TRGT_AP = case_when(
    TRGT_AP < 0.25 ~ "[0~0.25)",
    TRGT_AP < 0.5 ~ "[0.25~0.5)",
    TRGT_AP < 0.75 ~ "[0.5~0.75)",
    TRUE ~ '[0.75~1]')) %>%
  count(check,TRGT_AP) %>% 
  pivot_wider(names_from = TRGT_AP,values_from = n) %>%
  rename(match = check) %>%
  select(match,"[0~0.25)","[0.25~0.5)","[0.5~0.75)",'[0.75~1]') -> d

head(d)
addWorksheet(wb, "simpleSTR_AP")
writeData(wb, "simpleSTR_AP", d)



head(basic.count.STRallele)
basic.count.STRallele %>%
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
  count(check,STR_length,TRGT_AP) %>%
  pivot_wider(names_from = STR_length,values_from = n) %>%
  rename(match = check) %>% mutate(TRGT_AP = factor(TRGT_AP,levels = c("[0~0.25)", "[0.25~0.5)", "[0.5~0.75)", "[0.75~1]"))) %>%
  rename(AP = TRGT_AP) %>%
  arrange(match, AP) %>%
  select(match,AP,"[0~50)","[50~100)","[100~150)",'[150~Inf)') -> e


addWorksheet(wb, "simpleSTR_STRlengthAP")
writeData(wb, "simpleSTR_STRlengthAP", e)

saveWorkbook(wb, "~/Desktop/KU/@research/STR/table/basic_STRlength_count_byGC_AP.xlsx", overwrite = TRUE)
