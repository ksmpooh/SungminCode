library(tidyverse)
setwd("/Users/ksmpooh/Desktop/KCDC/long_read/2022/VCF/marker_filter/rm_multi.allele")

short_raw <- read.table("HLA.shortread.GATK.raw_onlySNP_INFO.txt") %>% mutate(type = "short_raw") %>% mutate(Seq="short")
long_raw <- read.table("HLA.longread.GATK.raw_onlySNP_INFO.txt") %>% mutate(type = "long_raw") %>% mutate(Seq="long")

short_VQSR <- read.table("HLA.Shortread.Seq.GATK.recal.pass.onlySNP_SNPIDupdate_INFO.txt") %>% mutate(type = "short_VQSR") %>% mutate(Seq="short")
long_VQSR <- read.table("HLA.Longread.Seq.GATK.recal.pass.onlySNP_SNPIDupdate_checksampleID_INFO.txt") %>% mutate(type = "long_VQSR") %>% mutate(Seq="long")

#short_CNN <- read.table("HLA.shortread.1DCNNfilter_annotated_onlySNP_INFO.txt") %>% mutate(type = "short_CNN") %>% mutate(Seq="short")
#long_CNN <- read.table("HLA.longread.1DCNNfilter_annotated_onlySNP_INFO.txt") %>% mutate(type = "long_CNN")%>% mutate(Seq="long")
head(long_CNN)
head(short_raw)

#df <- short_raw

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
colnames(df)
library(geomtextpath)

df %>% #replace_na(list("0")) %>%
  filter(Seq == "short") %>%
  ggplot(aes(x=QD,fill =type)) + 
  geom_vline(xintercept=15, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=20, y=20000, label="Some text") + 
  #geom_textvline(label = "the strong cars", xintercept = 30, vjust = 1,angle=45) +
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
  theme(legend.position = "none")+
  geom_histogram(position="dodge")


colnames(df)
p1<-df %>% #replace_na(list("0")) %>%
    filter(Seq == "short") %>%
  ggplot(aes(x=QD,fill =type)) + 
  theme(legend.position = "none")+
  geom_vline(xintercept=2, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=7, y=20000, label="QD<2.0") + 
  geom_histogram(position="dodge")
p1
p2<-df %>% #replace_na(list("0")) %>%
    filter(Seq == "short") %>%
  ggplot(aes(x=QUAL,fill =type)) + 
  theme(legend.position = "none")+
  #xlim(0, 2000000) +
  geom_vline(xintercept=30, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=2000000, y=100000, label="QUAL<30.0") + 
  geom_histogram(position="dodge")
p2
p3<-df %>% #replace_na(list("0")) %>%
    filter(Seq == "short") %>%
  ggplot(aes(x=SOR,fill =type)) + 
  theme(legend.position = "none")+
  geom_vline(xintercept=3, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=6, y=50000, label="SOR>3.0") + 
  geom_histogram(position="dodge")
p3

p4<-df %>% #replace_na(list("0")) %>%
    filter(Seq == "short") %>%
  ggplot(aes(x=FS,fill =type)) + 
  theme(legend.position = "none")+
  geom_vline(xintercept=60, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=100, y=75000, label="FS>60.0") + 
  geom_histogram(position="dodge")
p4

p5<-df %>% #replace_na(list("0")) %>%
    filter(Seq == "short") %>%
  ggplot(aes(x=MQ,fill =type)) + 
  theme(legend.position = "none")+
  geom_vline(xintercept=40, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=50, y=40000, label="MQ<40.0") + 
  geom_histogram(position="dodge")

p6<-df %>% #replace_na(list("0")) %>%
    filter(Seq == "short") %>%
  ggplot(aes(x=MQRankSum,fill =type)) + 
  theme(legend.position = "none")+
  geom_vline(xintercept=-12.5, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=-20, y=30000, label="MQRankSum < -12.5") +
  geom_histogram(position="dodge")

p7<-df %>% #replace_na(list("0")) %>%
    filter(Seq == "short") %>%
  ggplot(aes(x=ReadPosRankSum,fill =type)) + 
  theme(legend.position = "none")+
  geom_vline(xintercept=-8, linetype = 'dotted', color='red', size = 1) + 
  annotate("text", x=-10, y=30000, label="ReadPosRankSum < -4.0") + 
  geom_histogram(position="dodge")

legend <- cowplot::get_legend(df %>% #replace_na(list("0")) %>%
                                filter(Seq == "short") %>%
                                ggplot(aes(x=QD,fill =type)) + 
                                geom_histogram(position="dodge")+
                                #theme(legend.title = element_blank())
                                scale_fill_discrete(name = "SNP"))

grid.newpage()
grid.draw(legend)

library(grid)
library(gridExtra)
library(ggplot2)

grid.arrange(p1, p2, p3,p4,p5,p6,p7,legend, ncol=4, nrow = 2, 
             top = textGrob("GATK pipeline SNP Info. distribution for Hardfiltering (short-read)"))
             #layout_matrix = rbind(c(1,2,3), c(4,5 ,6), c(7,7,7)),
             #widths = c(2.7, 2.7,2.7), heights = c(2.5, 2.5,0.3))

