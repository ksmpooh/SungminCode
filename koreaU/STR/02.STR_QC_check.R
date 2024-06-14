library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/")

revio <- read_table("Revio.STR.pbmm2_hg38_withunmapped_trgt_genotype_gangstr.sorted.merged_withall.mkdf.txt")
head(revio)
short <- read_table("Shortread.STR.goast_hg38_withunmapped_EH_genotype_gangstr.merged_withall.mkdf.txt")
head(short)
revio %>% filter(ALT == ".")

short %>% filter(FILTER == "PASS") %>% filter(ALT != ".") 
  
  
## trgt quality check
setwd("/Users/ksmpooh/Desktop/KU/@research/STR/trgt/Quality_check/")
flist = grep(list.files("./"),pattern = "txt", value=TRUE)
df <- NULL
for (i in flist) {
  tmp <- read_table(i)
  df <- rbind(df,tmp)
}
df <- rbind(df,tmp)
head(df)

df %>% select()
df %>% count(TRID) %>% dim
df %>% filter(AP == "0")
df %>% filter(AP %in% c(".","0",NA)) %>% count(AP,MC)
df %>% filter(is.na(MC))
df %>% filter(MC != "0") %>% mutate(AP = as.numeric(AP)) -> new_df


head(new_df)
new_df %>% select(ID,TRID,Allele,MC,AP) %>% mutate(AP_range = ifelse(AP > 0.9,"AP > 0.9",ifelse()))
new_df %>% group_by(ID) %>% summarise(mean_AP = mean(AP),mean_MC = mean(MC)) %>% #head()
  ggplot(aes(x=ID,y=mean_AP)) +
  geom_point()


new_df %>% select(ID,TRID,Allele,MC,AP) %>% 
  group_by(TRID) %>% summarise(mean_AP = mean(AP)) %>% head()
  ggplot(aes(x=TRID,y=mean_AP)) +
  geom_point() + 
  theme(axis.text.x = element_blank())
  
