#### STR bed vs JSON
library(tidyverse)
library(ggvenn)


library(jsonlite)
library(httr)

setwd("/Users/ksmpooh/Desktop/KU/@research/STR")


#df <-read.table("Common_trgt_repeat_add_eh.v5_w_gangstr.v13.polymorphic.txt")
#head(df)

eh <- read.table("eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.bed")

#nrow(eh)
#nrow(eh_json)
trgt <- rbind(read.table("pathogenic_repeats.hg38.bed"),read.table("repeat_catalog.hg38.bed"))
trgt_path <- read.table("pathogenic_repeats.hg38.bed")
trgt_ori <- read.table("repeat_catalog.hg38.bed")
head(trgt)

eh_json <- fromJSON("EH_hg38_pipeline/resources/eh.v5_w_gangstr.v13.polymorphic.json")
head(eh_json)
eh$theme <- "eh"
trgt$theme <- "trgt"

head(eh)
head(trgt)

eh %>% mutate(ID = str_split_fixed(V4,";",3)[,1]) %>% #head()
  mutate(ID = gsub("ID=","",ID)) %>% 
  mutate(ID_pattern = ifelse(grepl("chr",V4),"chr","ID")) %>% 
  mutate(posID = paste0(V1,"_",V2,"_",V3))-> eh

trgt %>% mutate(ID = str_split_fixed(V4,";",3)[,1]) %>% #head()
  mutate(ID = gsub("ID=","",ID)) %>%
  mutate(ID_pattern = ifelse(grepl("chr",V4),"chr","ID")) %>% 
  mutate(posID = paste0(V1,"_",V2,"_",V3)) -> trgt

eh %>% rbind(trgt) %>% mutate(ID = str_split_fixed(V4,";",3)[,1]) %>% #head()
  mutate(ID = gsub("ID=","",ID)) -> df


eh %>% %>% count(ID_pattern)
head(eh)

head(df)
colnames(df)
df$ID

a <- eh %>% filter(ID_pattern == "chr") %>% select(ID)
b <- trgt %>% filter(ID_pattern == "chr") %>% select(ID)
v <- list("EH" = a$ID, "TRGT" = b$ID)
ggvenn(v,text_size = 7) + ggtitle("Common Repeat Pattern (chr_start_end)")



a <- eh %>% filter(ID_pattern != "chr") %>% select(ID)
b <- trgt %>% filter(ID_pattern != "chr") %>% select(ID)
v <- list("EH" = a$ID, "TRGT" = b$ID)
ggvenn(v,text_size = 7) + ggtitle("Pathogenic Repeat Pattern (ID)")


a <- eh %>% select(ID)
b <- trgt  %>% select(ID)
v <- list("EH" = a$ID, "TRGT" = b$ID)
ggvenn(v,text_size = 7) + ggtitle("Repeat Pattern (Overall)")


eh %>% filter(ID_pattern == "ID") %>% filter(!(ID %in% trgt$ID)) %>% filter(posID %in% trgt$posID)

trgt %>% filter(ID_pattern == "ID") %>% filter(!(ID %in% eh$ID)) %>% filter(posID %in% eh$posID)
trgt %>% filter(posID %in% c("chr1_149390802_149390841","chr13_70139353_70139428"))



trgt %>% filter(ID_pattern == "ID") %>% filter((ID %in% eh$ID)) -> trgt_ID_uniq

merge(trgt_ID_uniq,eh, all.x = T,by="posID") %>% select(posID,ID.x,ID.y,V4.x,V4.y) -> a

head(eh)
eh %>% select(posID) %>% duplicated() -> dup_index
eh %>% filter(dup_index) -> dup_posID
eh %>% filter(posID %in% dup_posID$posID) %>% filter(ID_pattern == "chr") -> eh_dup_chr
eh %>% filter(posID %in% dup_posID$posID) %>% filter(ID_pattern != "chr") -> eh_dup_ID
dup_posID
merge(eh_dup_chr,eh_dup_ID,all = T,by="posID") %>% select(V4.x,V4.y) %>% head()
eh_dup_ID %>% select(ID,posID) -> a
eh %>% dim()
eh %>% filter(!(ID %in% eh_dup_ID$ID)) %>% dim()
eh %>% filter(!(ID %in% eh_dup_ID$ID)) -> eh_rmdup


#write.table(eh_rmdup, "eh.v5_w_gangstr.v13.polymorphic.JSONtoBED.rmdup.bed",col.names = F,row.names = F,quote = F,sep = "\t")


a <- eh_rmdup %>% filter(ID_pattern == "chr") %>% select(ID)
b <- trgt %>% filter(ID_pattern == "chr") %>% select(ID)
v <- list("EH" = a$ID, "TRGT" = b$ID)
ggvenn(v,text_size = 7) + ggtitle("Common Repeat Pattern (chr_start_end)")



a <- eh_rmdup %>% filter(ID_pattern != "chr") %>% select(ID)
b <- trgt %>% filter(ID_pattern != "chr") %>% select(ID)
v <- list("EH" = a$ID, "TRGT" = b$ID)
ggvenn(v,text_size = 7) + ggtitle("Pathogenic Repeat Pattern (ID)")


a <- eh_rmdup %>% select(ID)
b <- trgt  %>% select(ID)
v <- list("EH" = a$ID, "TRGT" = b$ID)
ggvenn(v,text_size = 7) + ggtitle("Repeat Pattern (Overall)")

eh_rmdup %>% rbind(trgt) %>% count(theme,ID_pattern) -> a

trgt %>% filter(ID_pattern == "ID") %>% count()
trgt %>% filter(ID_pattern == "ID") %>% filter(ID %in% eh_rmdup$ID) %>% count()
trgt %>% filter(ID_pattern == "ID") %>% filter(posID %in% eh_rmdup$posID) %>% count()

trgt %>% filter(ID_pattern == "ID") %>% filter(!(ID %in% eh_rmdup$ID)) %>% count()
trgt %>% filter(ID_pattern == "ID") %>% filter(!(ID %in% eh_rmdup$ID)) %>% filter(posID %in% eh_rmdup$posID)
trgt %>% filter(ID_pattern == "ID") %>% filter(!(ID %in% eh_rmdup$ID)) %>% filter(posID %in% eh_rmdup$posID) %>% select(ID)
