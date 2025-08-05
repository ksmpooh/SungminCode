#### to kim hyung woo 20250117
library(tidyverse)


setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_ori")
df <- read_table("AR.20250117.withMono.genotype.txt")
df_noMono <- read_table("AR.20250116.txt")
nonMono <- read.table("non_mono.snp.txt")
head(df)
ncol(df)
ncol(df_noMono)
head(df_noMono)
dim(df)
ref <-read_table("AR.20250117.withMono.genotype.snpinfo.txt",col_names = F)

head(ref)
dim(ref)

df_noMono_header <- colnames(df_noMono)
df_header <- colnames(df)
df_noMono_header
head(ref)
ref %>% mutate(header = paste0(X3,"_chr",str_replace_all(X1,":","\\."))) -> new_header_pre
#ref %>% mutate(header = paste0(X3,"_",X1)) -> new_header_pre

new_header <- c("KBA_ID_KR", "KBA_ID_KD")

# 루프를 통해 새로운 헤더 생성
for (i in new_header_pre$header) {
  new_header <- c(new_header, paste0(i, "_KR"), paste0(i, "_KD"))
}
colnames(df)
new_header <- intersect(new_header, colnames(df))
#head(new_header)
head(new_header)
df %>% select(new_header)
head(ref)
dim(ref)

head(nonMono)
ref %>% mutate(monomorphic = ifelse(X1 %in% nonMono$V1,0,1)) -> out
colnames(out) <- c("ID","Annotation","Gene","monomorphic_SNP")
head(out)
head(df)

ref <- read_table("~/Desktop/KCDC/transplantation/00.sampleInfo/JG.IDupdate.NIHtobCODE.txt")
head(df)
head(ref)
df %>% select(new_header) %>%
  left_join(ref %>% rename(KBA_ID_KR = KBA_ID,bCODE_R = bCODE)) %>%  
  left_join(ref %>% rename(KBA_ID_KD = KBA_ID,bCODE_D = bCODE)) %>% 
  select(ncol(df)+1,ncol(df)+2,3:ncol(df)) -> df_out

writexl::write_xlsx(df_out,"AR.20250117.xlsx")
writexl::write_xlsx(out,"AR.20250117.snpinfo.xlsx")


### score

df <- read.table("KR.KD.AR_2025.alleleMatching01.Score_Sum.txt",header = T)
head(df)
head(ref)
df %>% #select(new_header) %>%
  left_join(ref %>% rename(KBA_ID.KR = KBA_ID,bCODE_R = bCODE)) %>%  
  left_join(ref %>% rename(KBA_ID.KD = KBA_ID,bCODE_D = bCODE)) %>% 
  select(ncol(df)+1,ncol(df)+2,3:ncol(df)) -> df_out

head(df_out)
writexl::write_xlsx(df_out,"AR.20250117.missmatch.ScoreSum.xlsx")


## check

geno_df <- readxl::read_xlsx("AR.20250117.xlsx")
head(geno_df)

geno_df %>% left_join(df_out) %>% 
  select(matches("LIMS")) %>% mutate(a = ifelse(LIMS1_chr2.109276221.G.T_KR == LIMS1_chr2.109276221.G.T_KD,0,1),
                                     b = ifelse(LIMS1_chr2.109292461.A.G_KR == LIMS1_chr2.109292461.A.G_KR,0,1)) %>% filter(a+b == LIMS1_Score_SUM) -> a









#### to YS prof.yang 20250131


library(tidyverse)


setwd("/Users/ksmpooh/Desktop/KCDC/transplantation/allogenomic/2025_AR/AR_YS/")
df_YS <- read_table("KR.KD.YS_31samples.NKcell.genotype.txt")
df_YS_mono <- read_table("KR.KD.YS_31samples.NKcell.genotype.withmono.txt")
head(df_YS)
ncol(df_YS)
head(df_YS_mono)
ncol(df_YS_mono)
monolist <- read.table("monosnp.txt")
head(monolist)
ref <-read_table("../AR_ori/AR.20250117.withMono.genotype.snpinfo.txt",col_names = F)
head(ref)
dim(ref)

df_YS_header <- colnames(df_YS)
df_YS_mono_header <- colnames(df_YS_mono)
df_YS_mono_header

head(ref)
ref %>% mutate(header = paste0(X3,"_chr",str_replace_all(X1,":","\\."))) -> new_header_pre
#ref %>% mutate(header = paste0(X3,"_",X1)) -> new_header_pre
#new_header_pre
new_header <- c("KBA_ID_KR", "KBA_ID_KD")

# 루프를 통해 새로운 헤더 생성
for (i in new_header_pre$header) {
  new_header <- c(new_header, paste0(i, "_KR"), paste0(i, "_KD"))
}
colnames(df_YS_mono)
new_header
new_header <- intersect(new_header, colnames(df_YS_mono))
#head(new_header)
head(new_header)
df_YS_mono %>% select(new_header)
head(ref)
dim(ref)

table(df_YS_mono_header %in% df_YS_header)
head(df_YS)
colnames(df_Y)
head(ref)
head(monolist)

ref %>% mutate(monomorphic = ifelse(X1 %in% monolist$V1,1,0)) -> out
colnames(out) <- c("ID","Annotation","Gene","monomorphic_SNP")




writexl::write_xlsx(df_YS_mono %>% select(new_header),"KR.KD.YS_31samples.NKcell.genotype.20250131.xlsx")
writexl::write_xlsx(out,"KR.KD.YS_31samples.NKcell.genotype.20250131.snpinfo.xlsx")



ncol(df_YS)
ncol(df)
