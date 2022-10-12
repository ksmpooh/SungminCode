library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KCDC/long_read/2022/VCF/marker_filter/rm_multi.allele")

short_raw <- read.table("HLA.shortread.GATK.raw_onlySNP_INFO.txt") %>% mutate(type = "short_raw") %>% mutate(Seq="short")
long_raw <- read.table("HLA.longread.GATK.raw_onlySNP_INFO.txt") %>% mutate(type = "long_raw") %>% mutate(Seq="long")

short_VQSR <- read.table("HLA.Shortread.Seq.GATK.recal.pass.onlySNP_SNPIDupdate_INFO.txt") %>% mutate(type = "short_VQSR") %>% mutate(Seq="short")
long_VQSR <- read.table("HLA.Longread.Seq.GATK.recal.pass.onlySNP_SNPIDupdate_checksampleID_INFO.txt") %>% mutate(type = "long_VQSR") %>% mutate(Seq="long")

#short_CNN <- read.table("HLA.shortread.1DCNNfilter_annotated_onlySNP_INFO.txt") %>% mutate(type = "short_CNN") %>% mutate(Seq="short")
#long_CNN <- read.table("HLA.longread.1DCNNfilter_annotated_onlySNP_INFO.txt") %>% mutate(type = "long_CNN")%>% mutate(Seq="long")
head(long_CNN)


df <- rbind(short_raw,long_raw) %>% rbind(short_CNN) %>% rbind(long_CNN) %>% rbind(long_VQSR) %>% rbind(short_VQSR)
df <- rbind(short_raw,long_raw) %>% rbind(long_VQSR) %>% rbind(short_VQSR)
table(df$type)

#'%POS\t%QD\t%QUAL\t%SOR\t%FS\t%MQ\t%MQRankSum\t%ReadPosRankSum\n'
colnames(df) <- c('POS','QD','QUAL','SOR','FS','MQ','MQRankSum','ReadPosRankSum','type',"Seq")
colnames(df)
head(df)
str(df)

df %>% filter(is.na('MQRankSum'))
df %>% #replace_na(list("0")) %>%
  mutate_at(c('QD','MQRankSum','ReadPosRankSum'),as.numeric) %>%
  group_by(type) %>%
  pivot_longer(cols = 2:8,names_to = "QC",values_to = "Value") %>%
  #gather(key=type,value=Value) %>% head()
  #group_by(type) %>%
  ggplot(aes(x=Value, fill = type)) + 
  geom_histogram(position="dodge") + 
  facet_wrap(~QC)


df <- df %>% #replace_na(list("0")) %>%
  mutate_at(c('QD','MQRankSum','ReadPosRankSum'),as.numeric)

df %>% #replace_na(list("0")) %>%
  filter(Seq == "short") %>%
  ggplot(aes(x=QD,fill =type)) + 
  geom_histogram(position="dodge")
  
df %>% #replace_na(list("0")) %>%
  filter(Seq == "long") %>%
  ggplot(aes(x=QD,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
#  filter(Seq == "short") %>%
  ggplot(aes(x=QD,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
  filter(Seq == "long") %>%
  ggplot(aes(x=QD,fill =type)) + 
  geom_histogram(position="dodge")


colnames(df)
df %>% #replace_na(list("0")) %>%
  #  filter(Seq == "short") %>%
  ggplot(aes(x=QD,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
  #  filter(Seq == "short") %>%
  ggplot(aes(x=QUAL,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
  #  filter(Seq == "short") %>%
  ggplot(aes(x=SOR,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
  #  filter(Seq == "short") %>%
  ggplot(aes(x=FS,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
  #  filter(Seq == "short") %>%
  ggplot(aes(x=MQ,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
  #  filter(Seq == "short") %>%
  ggplot(aes(x=MQRankSum,fill =type)) + 
  geom_histogram(position="dodge")

df %>% #replace_na(list("0")) %>%
  #  filter(Seq == "short") %>%
  ggplot(aes(x=ReadPosRankSum,fill =type)) + 
  geom_histogram(position="dodge")

