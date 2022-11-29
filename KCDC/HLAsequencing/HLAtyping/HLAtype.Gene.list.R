library(tidyverse)
library(stringr)
setwd("~/")
#setwd("/Volumes/DATA/HLAreferencePanel/")


#df <-read.table("HLA.nomenclean.chped")
#df <-read.table("HLA.type.result.8genes.merged.4digit_forMAKEreference.txt",header = T)
#head(df)

#writexl::write_xlsx(df,"~/Desktop/KCDC/HLA_seq/HLAtypeinfo.xlsx")

df <- readxl::read_excel("~/Desktop/KCDC/HLA_seq/HLAtypeinfo.xlsx",sheet = 2,na = ":")
head(df)

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  mutate('value' = str_split_fixed(value,"\\*",2)[,2]) %>% 
  group_by(HLAgene) %>% #head()
  #filter(HLAgene == "HLA-DRB1") %>% #head()
  count(value) %>% arrange(-n) -> df
  

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-A") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pA
  #facet_wrap(~HLAgene,ncol = 1) -> p
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#table(df$)
df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-B") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-C") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pC

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DRB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDRB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DQA1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDQA

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DQB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDQB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DPA1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDPA

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DPB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDPB

pA
pB
pC
pDRB
pDQA
pDQB
pDPA
pDPB

library(grid)
grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB, ncol=1) -> p
             #top = textGrob("GATK pipeline SNP Info. distribution for Hardfiltering (short-read)"))
png('~/Desktop/test.png',height = 5000)
print(p)
dev.off()

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-A") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  coord_flip()  #-> pDPB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  mutate('value' = str_split_fixed(value,"\\*",2)[,2]) %>% 
  group_by(HLAgene) %>% #head()
  #filter(HLAgene == "HLA-DRB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #facet_wrap(~HLAgene,ncol=2)



df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>%count(value)
  group_by(HLAgene,value) %>% count(HLAgene)
  summarise(n=n()) %>% arrange(HLAgene,-n)
  
#  #writexl::write_xlsx("~/Desktop/KCDC/HLA_seq/HLAtype.freq.xlsx")
