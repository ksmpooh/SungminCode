### 8K rare compore 130K NC paper rare sup4

library(tidyverse)
library(stringr)

setwd("/Users/ksmpooh/Desktop/KCDC/imputation/8K/asoo")
setwd("/ADATA/smkim/KBA_130K/06.asso/merge")
ref <- readxl::read_xlsx("~/Downloads/41467_2022_34163_MOESM3_ESM.xlsx",skip = 3,sheet = 8)
ref <- readxl::read_xlsx("/ADATA/smkim/KBA_130K/41467_2022_34163_MOESM3_ESM.xlsx",skip = 3,sheet = 8)
head(ref)

phenos = c('HDL_z','LDL_z','TC_z','TG_logz')

df <- read_table("HDL_z_asso.allCHR.txt")
colnames(ref)

ref %>% select(TRAIT,Novelty,rsID,CHR,BP,Effect...6,Other,MAF...10,`PVALUE...14`,Effect...12) %>%
  filter(TRAIT %in% c("HDL","LDL","TC","TG")) -> ref
head(ref)
ref$chr_pos <- paste0(as.character(ref$CHR),":", as.character(ref$BP))
df %>% mutate(chr_pos = paste0(str_split_fixed(SNP,":",4)[,1],":",str_split_fixed(SNP,":",4)[,2])) %>% filter(MAF < 0.01) -> df_rare
colnames(df_rare) <- c("new_beta","new_PVALUE","new_SNP","new_MAF","chr_pos")





phenos = c('HDL_z','LDL_z','TC_z','TG_logz')
ref %>% select(TRAIT,Novelty,rsID,CHR,BP,Effect...6,Other,MAF...10,`PVALUE...14`,Effect...12) %>%
  filter(TRAIT %in% c("HDL","LDL","TC","TG")) -> ref
head(ref)
ref$chr_pos <- paste0(as.character(ref$CHR),":", as.character(ref$BP))




target_trait <- c("HDL","LDL","TC","TG")
count = 0
out <- NULL
for (pheno in phenos) {
  count = count + 1
  df <- read_table(paste0(pheno,"_asso.allCHR.txt"))
  df %>% mutate(chr_pos = paste0(str_split_fixed(SNP,":",4)[,1],":",str_split_fixed(SNP,":",4)[,2])) -> df
  colnames(df) <- c("new_beta","new_PVALUE","new_SNP","new_MAF","chr_pos")
  ref_tmp <- ref %>% filter(TRAIT == target_trait[count])
  out <- rbind(out,ref_tmp %>% left_join(df))
}

write.table(out,"NCpaper130K_vs130K.8Kref.lipid.asso.txt",col.names = T,row.names = F,quote = F,sep = "\t")


#### after DP
df <- read.table("NCpaper130K_vs130K.8Kref.lipid.asso.txt",header = T)
head(df)          
writexl::write_xlsx(df,"NCpaper130K_vs130K.8Kref.lipid.asso.xlsx")
              



