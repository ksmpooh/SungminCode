## Repeat unit 비교
library(tidyverse)

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/Quality_check")
####
tmp <- read.table("NIH20N2000078.EH.qcmt.txt")
head(tmp)

tmp %>% select(-V1,-V2,-V9,-V10,-V11,-V12) %>% filter(V7 == 'PASS') %>%  select(-V7) %>%
  mutate(allele1=str_split_fixed(V8,"/",2)[,1],allele2=str_split_fixed(V8,"/",2)[,2]) %>%
  mutate(GT1 = str_split_fixed(V6,",",2)[,1],GT2 = str_split_fixed(V5,",",2)[,2]) %>% select(-V5,-V7) %>% head()

head(tmp)
tmp %>% select(3:11) -> tmp
tmp %>% select(V3,V4,V6) %>% filter(V6 == 'PASS') %>%  select(-V6) %>% count(V4)
####


flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
#  tmp %>% filter(X6 == "PASS") %>% select(X3,X4,X5) -> tmp
  tmp %>% select(3:12) -> tmp
  tmp$ID <- str_replace(i,".EH.qcmt.txt","")

  df <- rbind(df,tmp)
}
head(df)
#colnames(df) <- c("chr","pos","STR_ID","RU","ALT","REF","FILTER","GT","type","ADSP","ADFL","ADIR","ID")
#colnames(df) <- c("STR_ID","RU","ALT","ID")
eh <- df
head(eh)

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/Quality_check")
####
tmp <- read.table("NIH23F1013274.pbmm2_hg38_withunmapped_trgt_genotype_gangstr.sorted.MC_AP_AM.txt",header = T)
head(tmp)
tmp %>% filter(AP == 1)  %>% head()
tmp %>% filter(TRID == "chr1_151316_151328")
tmp %>% filter(AP == 1) %>% select(ID,TRID,MOTIFS,Allele,MC) %>% pivot_wider(names_from = Allele,values_from = MC) %>% #head
  filter(is.na(allele_1))
tmp %>% select(1,4,5,6,8,9,10,11) -> tmp

###
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
count = 0
for (i in flist) {
  count = count + 1
  print(count)
  tmp <- read_table(i)
#  tmp <- fread(i)
  tmp$ID <- str_replace(str_replace(i,".pbmm2_hg38_withunmapped_trgt_genotype_gangstr.sorted.MC_AP_AM.txt",""),"_sorted_trgt_genotype_gangstr.sorted.MC_AP_AM.txt","")
  tmp %>% select(1,4,5,6,8,9,10,11) -> tmp
  df <- rbind(df,tmp)
}
head(df)

tgrt <- df
df <- NULL
head(tgrt)
head(eh)
colnames(tgrt)
colnames(eh) <- c("STR_ID","RU","REF","ALT","FILTER","GT","type","ADSP","ADFL","ADIR","ID")
#head(eh)
#save(eh,file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.info.RData")
#save(tgrt,file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.info.RData")

###
load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.info.RData")
load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.info.RData")


trgt <- tgrt
tgrt <- NULL

head(eh)
head(trgt)
trgt %>% select(ID,TRID,Allele,AP) #%>% filter(AP != ".") %>% #head()
  #mutate(AP = as.numeric(AP)) %>%
  #pivot_wider(names_from = Allele,values_from = AP) %>% mutate(AP_mean = mean(allele_1,allele_2))

trgt %>% pivot_wider(names_from = Allele,values_from = AP)
eh %>% filter()

trgt %>% filter(!(AP %in% c(".","NA"))) %>% filter(MC != 0) %>% select(-AM) %>% filter(AP >= 0.8) %>% #dim() 41913285
  count(TRID) -> trgt_str_allele_count

head(trgt_str_allele_count)
#colnames(eh)[1] <- "STR_ID"
#colnames(eh)[5] <- "FILTER"
colnames(eh)[5] <- "QC"
head(eh)
eh %>% filter(QC == "PASS") %>% count(STR_ID) -> eh_str_allele_count

head(trgt_str_allele_count)
head(eh_str_allele_count)


trgt_str_allele_count %>% filter(n == 132) %>% dim()
trgt_str_allele_count[trgt_str_allele_count$n == 132,]$TRID
eh_str_allele_count %>% filter(n == 66)
head(eh_str_allele_count)
eh_str_allele_count[eh_str_allele_count$n == 66,]$STR_ID

trgt_str_allele_count[trgt_str_allele_count$n == 132,]$TRID %in% eh_str_allele_count[eh_str_allele_count$n == 66,]$STR_ID
trgt_str_allele_count
#trgt_str_allele_count %>% filter(n == 132) %>% filter(TRID %in% eh_str_allele_count[eh_str_allele_count$n == 66,]$STR_ID) %>% write.table("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/TRID_common_EHpass_TRGTupper0.8.txt",col.names = T,row.names = F,quote = F,sep = "\t")

library(ggvenn)

x<- list(
  EH = eh_str_allele_count[eh_str_allele_count$n == 66,]$STR_ID,
  TRGT = trgt_str_allele_count[trgt_str_allele_count$n == 132,]$TRID
)

ggvenn(x,set_name_size = 5,
       text_size = 5
        )


##
#load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.info.RData")
#load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.info.RData")


#trgt <- tgrt
#tgrt <- NULL

common_str <- read.table("~/Desktop/KU/@research/STR/02.compare/STR_type/TRID_common_EHpass_TRGTupper0.8.txt",header = T)

head(common_str)
head(eh)
head(trgt)
#colnames(eh)[1] <- "STR"
eh %>% filter(STR_ID %in% common_str$TRID) %>% select(-FILTER) -> eh_pass_common
#save(eh_pass_common,file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.PASS.intersect.RData")

trgt %>% filter(TRID %in% common_str$TRID) -> trgt_pass_common
#save(trgt_pass_common,file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.upper0.8.intersect.RData")
  

## PASS
load(file="~/Desktop/KU/@research/STR/Rdata/shortread.eh.66sample.qc.PASS.intersect.RData")
load(file="~/Desktop/KU/@research/STR/Rdata/longread.tgrt.66sample.qc.upper0.8.intersect.RData")

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare")
ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)
head(ref)
ref %>% select(Revio,Illumina) -> ref

head(eh_pass_common,10) %>% 
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>%
  mutate(ALT = ifelse(ALT == ".",0,ALT))
  
  

head(eh_pass_common)
head(trgt_pass_common)
trgt_pass_common %>% filter(grepl(",",MOTIFS))

eh_pass_common %>% select(-type,-ADSP,-ADFL,-ADIR) %>% 
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>%
  mutate(ALT = ifelse(ALT == ".",0,ALT)) %>% mutate(GT = str_split_fixed(GT,"/",2),ALT= str_split_fixed(ALT,",",2)) -> eh_pass_common_gt


eh_pass_common_gt %>% mutate(STR1 = ifelse(GT[,1] == "1",ALT[,1],ifelse(GT[,1] == "2",ALT[,2],ifelse(GT[,1] == "0",REF,".")))) %>%
  mutate(STR2 = ifelse(GT[,2] == "1",ALT[,1],ifelse(GT[,2] == "2",ALT[,2],ifelse(GT[,2] == "0",REF,".")))) %>% 
  select(-ALT,-GT,-REF) %>% mutate(STR1 = as.integer(STR1),STR2 = as.integer(STR2)) %>%
  mutate(EU_STR1 = ifelse(STR1 < STR2, STR1,STR2),EU_STR2 = ifelse(STR1 < STR2, STR2,STR1)) %>% select(-STR1,-STR2)-> eh_pass_common_gt

#rm(eh_pass_common_gt1)
head(eh_pass_common_gt)
#eh_pass_common_gt %>% filter(GT[,1] == ".")
#eh_pass_common_gt %>% filter(GT[,2] == ".")
head(trgt_pass_common)

trgt_pass_common %>% select(-AP,-AM,-STRUC) %>% pivot_wider(names_from = Allele,values_from = MC) %>% 
  mutate(TRGT_STR1 = ifelse(allele_1 < allele_2, allele_1,allele_2),TRGT_STR2 = ifelse(allele_1 < allele_2, allele_2,allele_1)) %>% 
  select(-allele_1,-allele_2)-> trgt_pass_common_gt

eh_pass_common_gt %>% count(EU_STR1 <= EU_STR2)
trgt_pass_common_gt %>% count(TRGT_STR1 <= TRGT_STR2)

head(eh_pass_common_gt)
head(trgt_pass_common_gt)


colnames(eh_pass_common_gt) <- c("STR","MOTIFS","EH_ID","EH_STR1","EH_STR2")
colnames(trgt_pass_common_gt) <- c("ID","STR","MOTIFS","TRGT_STR1","TRGT_STR2")

colnames(ref) <- c("ID","EH_ID")

eh_pass_common_gt %>% left_join(ref) %>% select(-EH_ID) %>% left_join(trgt_pass_common_gt) -> df
head(df)
#이거 문제 있어서.... 다른걸로 체크
#write_tsv(df,"~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReeaptnumber.txt")

#df <- read_table("~/Desktop/KU/@research/STR/02.compare/STR.TRGT_0.8upper.EH_pass.common.merge.onlyReeaptnumber.txt")
head(df)
df %>% filter(ID == 'NIH23F1724125',STR=="chr1_144527_144575")
df %>% #group_by(ID,STR) %>% 
  mutate(STR1 = TRGT_STR1 - EH_STR1,STR2 = TRGT_STR2 - EH_STR2) %>%
  select(-TRGT_STR1,-TRGT_STR2,-EH_STR1,-EH_STR2) -> df_diff


head(df_diff)

df_diff %>% filter(STR1 !=0 | STR2 != 0) %>% pivot_longer(cols = c("STR1","STR2"),names_to = "Allele") %>%
  group_by(ID,STR) %>% head


df_diff %>% 



