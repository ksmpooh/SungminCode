library(tidyverse)
library(stringr)
library(readxl)
library(ggbreak)
library(ggpubr)
library(dplyr)
library(patchwork) # To display 2 charts together
library(hrbrthemes)


setwd("~/Desktop/KCDC/HLA_seq/")

df <- read_xlsx("variant.xlsx",sheet = 2)
head(df)

colnames(df)
df %>% filter(TOOL == "GATK") %>% #head() 
  #filter(Seqeuncing )
  select(Seqeuncing,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`,...12) %>%  
  mutate(QC = ifelse(is.na(...12),QC,paste0(QC,"_",...12))) %>% #head()
  group_by(Seqeuncing) %>%
  pivot_longer(cols = c(RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%#  head()
  ggplot(aes(x=type,y=Val,color = QC)) +
  geom_point() +
  geom_line() + 
  #geom_bar() +
  facet_grid(~Seqeuncing)
  

df %>% filter((is.na(...12))) %>%
  select(Seqeuncing,TOOL,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  pivot_longer(cols = c(RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%  #head()
  ggplot(aes(x=type,y=Val,color = QC)) +
  geom_point() +
  geom_line() + 
  facet_grid(TOOL~Seqeuncing)


#aes(x = key,fill=factor(value,levels = c('Caucasian','Asian','Slavic','Black'))))

df %>% filter((is.na(...12))) %>% #count(TOOL,QC)
  #select(Seqeuncing,TOOL,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  select(Seqeuncing,TOOL,QC,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  pivot_longer(cols = c(SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%  #head()
  mutate(Method = paste0(TOOL,"_",QC)) %>% #head()
  mutate(Method = str_replace(Method, "Deepvariant-GLnexus", "DV")) %>% #count(Method)
  #ggplot(aes(x=Method,y=Val,fill = type)) +
  ggplot(aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),
             y=Val,fill = factor(type,levels=c("Multiallelic SNP site","Multiallelic site","INDEL","SNP")))) +
  geom_bar(stat = "identity") + 
  theme(legend.title = element_blank(),legend.position="bottom") +
  xlab(element_blank()) + ylab(element_blank()) + 
  #scale_fill_discrete(guide="none") +
  facet_grid(~Seqeuncing) #->a1

  #geom_point() #+
  #geom_line() + 


df %>% filter((is.na(...12))) %>% #count(TOOL,QC)
  #select(Seqeuncing,TOOL,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  select(Seqeuncing,TOOL,QC,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  pivot_longer(cols = c(SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`),names_to = 'type',values_to = 'Val') %>%  #head()
  mutate(Method = paste0(TOOL,"_",QC)) %>% #head()
  mutate(Method = str_replace(Method, "Deepvariant-GLnexus", "DV")) %>% #count(Method)
  #ggplot(aes(x=Method,y=Val,fill = type)) +
  ggplot(aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),
             y=Val,fill = factor(type,levels=c("Multiallelic SNP site","Multiallelic site","INDEL","SNP")))) +
  geom_bar(position = 'fill',stat = "identity") + 
  theme(legend.title = element_blank(),legend.position="bottom") +
  xlab(element_blank()) + ylab(element_blank()) + 
  facet_grid(~Seqeuncing) ->a2


ggarrange(a1,a2,ncol = 1,nrow = 2,heights = c(8,9))
#################3
df %>% filter((is.na(...12))) %>% #count(TOOL,QC)
  #select(Seqeuncing,TOOL,QC,RECORDS,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`) %>%  #head()
  select(Seqeuncing,TOOL,QC,SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`,`TS/TV`,`TS/TV(1st ALT)`) %>%  #head()
  pivot_longer(cols = c(SNP,INDEL,`Multiallelic site`,`Multiallelic SNP site`,`TS/TV`,`TS/TV(1st ALT)`),names_to = 'type',values_to = 'Val') %>%  #head()
  mutate(Method = paste0(TOOL,"_",QC)) %>% #head()
  mutate(Method = str_replace(Method, "Deepvariant-GLnexus", "DV")) -> df1#%>% #count(Method)

head(df1)
head(df1)
str(df1)
df1 %>% filter(grepl("TS",type))
ggplot() + 
  geom_bar(data = df1 %>% filter(!grepl("TS",type)),
           aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),
               y=Val,fill = factor(type,levels=c("Multiallelic SNP site","Multiallelic site","INDEL","SNP"))),stat = "identity") + 
  geom_point(data = df1 %>% filter(grepl("TS",type)),aes(x=Method,y=Val*200000,color=type)) + 
  #scale_fill_manual(values=c("TS/TV"="#002955" ,"TS/TV(1st ALT)"="#074ca1")) +
  #geom_point(data = df1 %>% filter(grepl("TS",type)),aes(x=Method,y=Val*200000,color=c(rgb(0,0,0,0),rgb(0,0,1,0.5)))) + 
  scale_y_continuous(name = 'Count (Bar)',sec.axis = sec_axis(~ ./200000,name = 'TS/TV ratio (Point)'))+
  #scale_color_discrete_qualitative(palette = "Cold") + 
  #theme_ipsum() + 
  theme(legend.title = element_blank(),legend.position="bottom") +
  xlab(element_blank()) + ylab(element_text("# of variant"))+
  facet_grid(~Seqeuncing) #->a1



#ggplot(aes(x=Method,y=Val,fill = type)) +
  ggplot() +
  geom_bar(aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),
               y=Val,fill = factor(type,levels=c("Multiallelic SNP site","Multiallelic site","INDEL","SNP"))),stat = "identity") + 
  #geom_line(aes(x=filter(grepl("TS",type)),y=Val,color=filter(grepl("TS",type)))) +
    geom_line(aes(x=filter(grepl("TS",type)),y=Val,color=filter(grepl("TS",type)))) +
  theme(legend.title = element_blank(),legend.position="bottom") +
  xlab(element_blank()) + ylab(element_blank()) + 
  #scale_fill_discrete(guide="none") +
  facet_grid(~Seqeuncing) #->a1

#geom_point() #+
#geom_line() + 
