library(tidyverse)
library(writexl)
setwd("/Volumes/DATA/HLA_seq/annotation/")

a <- read.csv("long_VEP_summary_dataprocessing.txt")
b <- read.csv("short_VEP_summary_dataprocessing.txt")

head(a)

a_vari <- a %>% filter(grepl("Variant",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2]) 
b_vari <- b %>% filter(grepl("Variant",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2])


a_mt <- a %>% filter(grepl("most",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2]) 
b_mt<- b %>% filter(grepl("most",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2])

a_all<- a %>% filter(grepl("all",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2]) 
b_all<- b %>% filter(grepl("all",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2])

a_coding <- a %>% filter(grepl("Coding",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2]) 
b_coding <- b %>% filter(grepl("Coding",index)) %>% mutate(type = str_split_fixed(index,'-',2)[,2])
head(a)
head(b)

vari <- merge(a_vari[,c(3,2)],b_vari[,c(3,2)])
mt <- merge(a_mt[,c(3,2)],b_mt[,c(3,2)])
all <- merge(a_all[,c(3,2)],b_all[,c(3,2)])
coding <- merge(a_coding[,c(3,2)],b_coding[,c(3,2)])


write_xlsx(list("General" = vari,"most" = mt,"all" = all,"coding" = coding),"~/Desktop/KCDC/HLA_seq/HLAseq_annotation.xlsx",)
head(anno_out)
str(anno_out)
anno_out$noInHanREF_count <- as.numeric(anno_out$noInHanREF_count)
anno_out$ALL_count <- as.numeric(anno_out$ALL_count)
anno_out$ALL_per <- anno_out$ALL_count/sum(anno_out$ALL_count) * 100
anno_out$noInHanREF_per <- anno_out$noInHanREF_count/sum(anno_out$noInHanREF_count) * 100
anno_out <- anno_out[,c(1,2,4,3,5)]
#write.csv(anno_out,"anno_Result.csv")
anno_out %>% pivot_longer(ALL_count,)

ggplot(anno_out,aes(x=type,color=noInHanREF_count)) +
  geom_histogram()
