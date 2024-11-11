### chrX

library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/chrX/raw_chrX")
sex_info <- read_table("~/Desktop/KCDC/pangenome/00.datacheck/Revio.WGS.sex.info.txt")
head(sex_info)
## TGRT quality check
#system.time(tmp <- read_table(i))
#system.time(tmp <- fread(i))
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
venn_rawdata <- read.table("~/Desktop/KU/@research/STR/figure/STR_DB_withQC_PASSofEH_PASSofTRGT_forVenn_v2.txt",header = T)
venn_rawdata %>% na.omit() %>% select(STR_DB) -> common_STR

head(common_STR)
head(tmp)
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
count = 0
for (i in flist) {
  #tmp <- read_table(i)
  count = count + 1
  print(count)
  tmp <- fread(i)
  tmp$ID <- i
  df <- rbind(df,tmp)
}

#df <- rbind(df,tmp)
head(df)


df %>% filter(!str_detect(ID,"NIH23J3904558"))  %>% filter(TRID %in% common_STR$STR_DB) %>% #head()
  mutate(across(MC:AM, ~ ifelse(. == ".", -1, .))) -> df


df %>% filter(str_detect(ID,"normal")) %>% filter(str_detect(TRID,"chr")) -> df_normalSTR_chrX
#head(df_normalSTR_chrX)

df %>% filter(str_detect(ID,"patho")) %>% filter(str_detect(MOTIFS,",")) -> df_pathoSTR_complex_chrX
df %>% filter(str_detect(ID,"patho")) %>% filter(!str_detect(MOTIFS,",")) -> df_pathoSTR_simple_chrX

df_normalSTR_chrX %>% rbind(df_pathoSTR_simple_chrX) %>% mutate(ID = str_split_fixed(ID,"\\.",5)[,1]) %>% #select(ID)
  write.table("../chrX.trgt.simpleSTR.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")



df_pathoSTR_complex_chrX%>% mutate(ID = str_split_fixed(ID,"\\.",5)[,1]) %>% #select(ID)
  write.table("../chrX.trgt.complexSTR.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")

############

setwd("/Users/ksmpooh/Desktop/KU/@research/STR/eh/Quality_check_chrX/")

#sex_info %>% filter(Revio == "NIH23J3904558")
complex_str <- read.table("~/Desktop/KU/@research/STR/figure/extra_info/complex_STR.v2.txt",header = T)
head(complex_str)
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
flist
df <- NULL
for (i in flist) {
  tmp <- read_table(i,col_names = FALSE)
  #  tmp %>% filter(X6 == "PASS") %>% select(X3,X4,X5) -> tmp
  #tmp %>% select(3:12) -> tmp
  #tmp$ID <- str_replace(i,".EH.qcmt.txt","")
  tmp$ID <- i
  
  df <- rbind(df,tmp)
}
head(df)
colnames(df) <- c("CHROM","POS","STR_ID","RU","REF","ALT","FILTER","GT","type","ADSP","ADFL","ADIR","LC","filename")
#colnames(df) <- c("STR_ID","RU","ALT","ID")
head(df)
head(sex_info)

ref <- read.table("~/Desktop/KCDC/pangenome/00.datacheck/KBA.Long_Revio_Nanopore_short.IDmatchinagtable.txt",header = T)

ref %>% select(Revio,Illumina) -> ref
head(ref)
colnames(ref) <- c("ID","EH_ID")

df %>% mutate(EH_ID = str_split_fixed(filename,"\\.",4)[,1]) %>% filter(EH_ID != "NIH20N2594890") %>%
  filter(str_detect(CHROM,"chrX")) %>% filter(EH_ID %in% sex_info[sex_info$sex == "M",]$Illumina) %>% left_join(ref) %>% select(-EH_ID) -> df_pro_chrX
head(df_pro_chrX)


df_pro_chrX %>% filter(str_detect(filename,"patho")) -> df_pro_chrX_patho



df_pro_chrX_patho %>% filter(STR_ID %in% complex_str$STR_ID) -> df_pro_chrX_patho_complexSTR
df_pro_chrX_patho %>% filter(!(STR_ID %in% complex_str$STR_ID)) %>% filter(STR_ID %in% common_STR$STR_DB)-> df_pro_chrX_patho_simpleSTR
head(df_pro_chrX_patho_simpleSTR)

#df_pro_chrX %>% filter(!str_detect(filename,"patho")) %>% filter(STR_ID %in% common_STR$STR_DB) %>% dim()
#df_pro_chrX %>% filter(!str_detect(filename,"patho")) %>% filter(STR_ID %in% common_STR$STR_DB) %>% duplicated() %>% dim()
df_pro_chrX %>% filter(!str_detect(filename,"patho")) %>% filter(STR_ID %in% common_STR$STR_DB) %>% filter(str_detect(STR_ID,"chr")) -> df_pro_chrX_normal
head(df_pro_chrX_normal)
df_pro_chrX_normal %>% filter(!str_detect(STR_ID,"chr"))

df_pro_chrX_normal %>% rbind(df_pro_chrX_patho_simpleSTR) %>% select(-filename) %>% select(-ADSP,-ADFL,-ADIR) %>% 
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>% #filter(GT == 0) %>% count(ALT)#count(GT != '0' & ALT != ".")
  mutate(STR1 = ifelse(GT == 0,REF,ALT)) %>% select(-ALT,-FILTER,-GT) %>% 
  write.table("/Users/ksmpooh/Desktop/KU/@research/STR/eh/chrX/chrX.eh.simpleSTR.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")


df_pro_chrX_patho_complexSTR %>% select(-filename) %>% select(-ADSP,-ADFL,-ADIR) %>% 
  mutate(ALT = str_replace_all(ALT,"STR",""),ALT = str_replace_all(ALT,"<",""),ALT = str_replace_all(ALT,">","")) %>% #filter(GT == 0) %>% count(ALT)#count(GT != '0' & ALT != ".")
  mutate(STR1 = ifelse(GT == 0,REF,ALT)) %>% select(-ALT,-FILTER,-GT) %>% #head()
  write.table("/Users/ksmpooh/Desktop/KU/@research/STR/eh/chrX/chrX.eh.complexSTR.processing.txt",col.names = T,row.names = F,quote = F,sep = "\t")


### chrX patho ccomplex

eh_pro_chrX_patho_complexSTR <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/eh/chrX/chrX.eh.complexSTR.processing.txt")
head(eh_pro_chrX_patho_complexSTR)

trgt_chrX_patho_complexSTR <- read_table("~/Desktop/KU/@research/STR/trgt/chrX/chrX.trgt.complexSTR.processing.txt")
final_ref <- read.table("~/Desktop/KU/@research/STR/db/Final.ref.20240829.txt",header=T)
head(final_ref)
final_ref %>% filter(str_detect(MOTIFS,",")) %>% select(ID,MOTIFS) %>% separate_rows(MOTIFS, sep = ",") %>% 
  mutate(new_ID = paste(ID, MOTIFS, sep = "_")) -> complex_ID
colnames(complex_ID)[1] <- c("STR_ID")

eh_pro_chrX_patho_complexSTR %>% mutate(new_ID = ifelse(str_detect(STR_ID,"_"),STR_ID,paste(STR_ID, RU, sep = "_"))) %>%
  mutate(STR_ID = str_split_fixed(new_ID,"_",2)[,1]) -> eh_pro_chrX_patho_complexSTR_merge
head(eh_pro_chrX_patho_complexSTR_merge)
head(complex_ID)
head(eh_complex_pass)
head(trgt_complex_pass)
trgt_chrX_patho_complexSTR %>% select(ID,TRID,MOTIFS,Allele,MC)
trgt_chrX_patho_complexSTR %>% separate_rows(MOTIFS, sep = ",") %>% mutate(new_ID = paste0(TRID,"_",MOTIFS)) %>% select(-MC) -> trgt_chrX_patho_complexSTR_MOTIFS
trgt_chrX_patho_complexSTR %>% #separate_rows(MOTIFS, sep = ",") %>% mutate(new_ID = paste0(TRID,"_",MOTIFS)) %>% 
  separate_rows(MC, sep = "_") %>% select(MC) -> trgt_chrX_patho_complexSTR_MC
head(trgt_chrX_patho_complexSTR_MC)
#trgt_complex_pass %>% separate_rows(MC, sep = ",") %>% mutate(new_ID = paste0(TRID,"_",MOTIFS)) %>% select(-MC) -> trgt_complex_pass_MOTIFS

#select(-MOTIFS) -> trgt_complex_pass_MC
head(trgt_chrX_patho_complexSTR_MOTIFS)
head(trgt_chrX_patho_complexSTR_MC)
dim(trgt_chrX_patho_complexSTR_MOTIFS)
dim(trgt_chrX_patho_complexSTR_MC)

trgt_chrX_patho_complexSTR_MOTIFS %>% cbind(trgt_chrX_patho_complexSTR_MC) %>% rename(STR_ID = TRID) %>% select(-STRUC,-ALT) -> trgt_chrX_patho_complexSTR_MC_merge

head(trgt_chrX_patho_complexSTR_MC_merge)
head(eh_pro_chrX_patho_complexSTR_merge)
eh_pro_chrX_patho_complexSTR_merge %>% select(-RU,-REF) %>% rename(EH_STR = STR1) %>%
  left_join(trgt_chrX_patho_complexSTR_MC_merge %>% select(-CHROM,-POS,-Allele) %>% rename(TRGT_STR = MC)) %>% 
  write.table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/chrX/eh_trgt_complexSTR_prep_chrX.txt",col.names = T,row.names = F,quote = F,sep = "\t")

#eh_pro_chrX_patho_complexSTR_merge %>% count(ID %in% trgt_chrX_patho_complexSTR_MC_merge$ID)

#"/Users/ksmpooh/Desktop/KU/@research/STR/figure/chrX/"
### chrX simple

eh_chrX_simpleSTR <- read_table("/Users/ksmpooh/Desktop/KU/@research/STR/eh/chrX/chrX.eh.simpleSTR.processing.txt")
trgt_chrX_simpleSTR <- read_table("~/Desktop/KU/@research/STR/trgt/chrX/chrX.trgt.simpleSTR.processing.txt")
head(eh_chrX_simpleSTR)
head(trgt_chrX_simpleSTR)
table(duplicated(eh_chrX_simpleSTR))
table(duplicated(trgt_chrX_simpleSTR))
dim(eh_chrX_simpleSTR)
dim(trgt_chrX_simpleSTR)


trgt_chrX_simpleSTR %>% rename(STR_ID = TRID) %>% select(-CHROM,-POS,-Allele,-STRUC,-ALT) %>% rename(TRGT_STR = MC) %>% head()
eh_chrX_simpleSTR %>% select(-RU,-REF) %>% rename(EH_STR = STR1) %>% head()

eh_chrX_simpleSTR %>% select(-RU,-REF) %>% rename(EH_STR = STR1) %>% #head()
  left_join(trgt_chrX_simpleSTR %>% rename(STR_ID = TRID) %>% select(-CHROM,-POS,-Allele,-STRUC,-ALT) %>% rename(TRGT_STR = MC)) %>% #head()
  write.table("/Users/ksmpooh/Desktop/KU/@research/STR/figure/chrX/eh_trgt_simpleSTR_prep_chrX.txt",col.names = T,row.names = F,quote = F,sep = "\t")
