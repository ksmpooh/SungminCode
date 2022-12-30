library(tidyverse)
library(stringr)
setwd("~/")
setwd("/Volumes/DATA/HLAreferencePanel/")


#df <-read.table("HLA.nomenclean.chped")
#df <-read.table("HLA.type.result.8genes.merged.4digit_forMAKEreference.txt",header = T)
#head(df)

#writexl::write_xlsx(df,"~/Desktop/KCDC/HLA_seq/HLAtypeinfo.xlsx")

df <- readxl::read_excel("~/Desktop/KCDC/HLA_seq/HLAtypeinfo.xlsx",sheet = 2,na = ":")
head(df)

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  mutate('value' = str_split_fixed(value,"\\*",2)[,2]) %>% 
  group_by(HLAgene) %>%  #head()
  count(value) %>% #head()#arrange(-n) -> df
  count(HLAgene) -> dfcount#%>%
  

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-A") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
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
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-C") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pC

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DRB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDRB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DQA1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDQA

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DQB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDQB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DPA1") %>% #head()
  count(value) %>% arrange(-n) %>% #head()
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
  #ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())+
  coord_flip()  -> pDPA

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DPB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=reorder(value,-n))) +
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
library(gridExtra)
library(cowplot)

layout_matrix
#grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB, ncol=1) -> p
grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB, nrow=6,
             layout_matrix = rbind(c(1,1),c(2,2)))

grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB, ncol = 2, 
             layout_matrix = cbind(c(1,1), c(2,2),c(3,3),c(4,4),c(5,6),c(7,8)))

grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB, nrow = 2, 
             layout_matrix = cbind(c(1,1), c(2,2),c(3,3),c(4,4),c(5,6),c(7,8)))

grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB,   ncol = 3, 
             layout_matrix = cbind(c(1,1,1), c(2,2,2),c(3,3,3),c(4,4,4),c(5,5,6),c(7,8,8)))

grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB,   widths = c(1,1,1,1,1,1),#labels = "AUTO",
             layout_matrix = cbind(c(1,1,1,1), c(2,2,2,2),c(3,3,3,3),c(4,4,4,4),c(5,5,6,6),c(7,8,8,8)))
cbind(c(1,1,1), c(2,2,2),c(3,3,3),c(4,4,4),c(5,5,6),c(7,8,8))
rbind(c(1,1,1), c(2,2,2),c(3,3,3),c(4,4,4),c(5,5,6),c(7,8,8))

##########
df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-A") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())-> pA
#facet_wrap(~HLAgene,ncol = 1) -> p
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#table(df$)
df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-B") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) ->pB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-C") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank())-> pC

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DRB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) -> pDRB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DQA1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) -> pDQA

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DQB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) -> pDQB

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DPA1") %>% #head()
  count(value) %>% arrange(-n) %>% #head()
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  #ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) -> pDPA

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-DPB1") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,-n),y=n,fill=reorder(value,-n))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) -> pDPB

grid.arrange(pA, pB, pC,pDRB,pDQA,pDQB,pDPA,pDPB,   widths = c(1,1,2),#labels = "AUTO",
             layout_matrix = rbind(c(1,1,1), c(2,2,2),c(3,3,3),c(4,4,4),c(5,5,6),c(7,8,8)))


###########3
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
  facet_wrap(~HLAgene,ncol=2)



df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% #count(value)
  group_by(HLAgene,value) %>% #count(HLAgene) 
  summarise(n=n()) %>% arrange(HLAgene,-n)
  
#  #writexl::write_xlsx("~/Desktop/KCDC/HLA_seq/HLAtype.freq.xlsx")
head(df)

#han <- read.csv("~/Desktop/KCDC/HLAimputation/HAN.ref/Han.HLAtyping.Result.4digit.csv",stringsAsFactors='TRUE')
han <- read_csv("~/Desktop/KCDC/HLAimputation/HAN.ref/Han.HLAtyping.Result.4digit.csv",col_names=TRUE)
head(han)
colnames(han)
grepl(".1",colnames(han))
han1 <- han %>% select(grep("\\.1",colnames(han)))
han2 <- han %>% select(grep("\\.2",colnames(han)))
head(han1)
head(han2)
colnames(han2) <- colnames(han1)
head(df)
han <- han1 %>% rbind(han2)
head(han)
#str_replace(colnames(han),"\\.1","")
colnames(han)
han %>% mutate(HLAtype_A.1 = paste0("A*",str_sub(HLAtype_A.1,1,-3),":",str_sub(HLAtype_A.1,-2,-1))) %>%
  mutate(HLAtype_B.1 = paste0("B*",str_sub(HLAtype_B.1,1,-3),":",str_sub(HLAtype_B.1,-2,-1))) %>%
  mutate(HLAtype_C.1 = paste0("C*",str_sub(HLAtype_C.1,1,-3),":",str_sub(HLAtype_C.1,-2,-1))) %>%
  mutate(HLAtype_DRB1.1 = paste0("DRB1*",str_sub(HLAtype_DRB1.1,1,-3),":",str_sub(HLAtype_DRB1.1,-2,-1))) %>%
  mutate(HLAtype_DPA1.1 = paste0("DPA1*",str_sub(HLAtype_DPA1.1,1,-3),":",str_sub(HLAtype_DPA1.1,-2,-1))) %>%
  mutate(HLAtype_DPB1.1 = paste0("DPB1*",str_sub(HLAtype_DPB1.1,1,-3),":",str_sub(HLAtype_DPB1.1,-2,-1))) %>%
  mutate(HLAtype_DQA1.1 = paste0("DQA1*",str_sub(HLAtype_DQA1.1,1,-3),":",str_sub(HLAtype_DQA1.1,-2,-1))) %>%
  mutate(HLAtype_DQB1.1 = paste0("DQB1*",str_sub(HLAtype_DQB1.1,1,-3),":",str_sub(HLAtype_DQB1.1,-2,-1))) %>% 
  pivot_longer(cols = colnames(han),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "HLAtype_", replacement = "HLA-")) %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "\\.1", replacement = "")) %>% #head()
  mutate('value' = str_split_fixed(value,"\\*",2)[,2]) %>% 
  #group_by(HLAgene) %>%  #head()
  group_by(HLAgene) -> a
  count(value) %>% #head()#arrange(-n) -> df
  count(HLAgene) -> hancount #%>%

head(hancount)
head(dfcount)

colnames(hancount)[2] <- "Han Chineses" 
colnames(dfcount)[2] <- "Korean" 
'
geom_bar(data = df1 %>% filter(!grepl("TS",type)),
         aes(x=factor(Method,levels=c("DV_unfiltered","DV_filtered","GATK_unfiltered","GATK_VSQR","GATK_Hardfiltering")),
             y=Val,fill = factor(type,levels=c("Multiallelic SNP site","Multiallelic site","INDEL","SNP"))),stat = "identity") + 
            
' 
#pivot_longer(cols = colnames(han),names_to = "HLAgene") %>% na.omit() %>% #head()
dfcount
a <- dfcount %>% left_join(hancount)
  
dfcount %>% left_join(hancount) %>% #head()
  pivot_longer(cols = c("Korean","Han Chineses"),names_to = "Ethnity") %>% #head()
  ggplot(aes(x=HLAgene,y=value,fill=factor(Ethnity,levels = c("Korean","Han Chineses")))) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_discrete(name = "Ethnity") +
  theme(legend.position="bottom") +
  ylab(element_blank()) + xlab(element_blank())
  #geom_bar(stat="identity")
  
  
han %>% mutate(HLAtype_A.1 = paste0("A*",str_sub(HLAtype_A.1,1,-3),":",str_sub(HLAtype_A.1,-2,-1))) %>%
  mutate(HLAtype_B.1 = paste0("B*",str_sub(HLAtype_B.1,1,-3),":",str_sub(HLAtype_B.1,-2,-1))) %>%
  mutate(HLAtype_C.1 = paste0("C*",str_sub(HLAtype_C.1,1,-3),":",str_sub(HLAtype_C.1,-2,-1))) %>%
  mutate(HLAtype_DRB1.1 = paste0("DRB1*",str_sub(HLAtype_DRB1.1,1,-3),":",str_sub(HLAtype_DRB1.1,-2,-1))) %>%
  mutate(HLAtype_DPA1.1 = paste0("DPA1*",str_sub(HLAtype_DPA1.1,1,-3),":",str_sub(HLAtype_DPA1.1,-2,-1))) %>%
  mutate(HLAtype_DPB1.1 = paste0("DPB1*",str_sub(HLAtype_DPB1.1,1,-3),":",str_sub(HLAtype_DPB1.1,-2,-1))) %>%
  mutate(HLAtype_DQA1.1 = paste0("DQA1*",str_sub(HLAtype_DQA1.1,1,-3),":",str_sub(HLAtype_DQA1.1,-2,-1))) %>%
  mutate(HLAtype_DQB1.1 = paste0("DQB1*",str_sub(HLAtype_DQB1.1,1,-3),":",str_sub(HLAtype_DQB1.1,-2,-1))) %>% 
  pivot_longer(cols = colnames(han),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "HLAtype_", replacement = "HLA-")) %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "\\.1", replacement = "")) %>% #head()
  mutate('value' = str_split_fixed(value,"\\*",2)[,2]) %>% 
  mutate(Ethnity = "Han Chineses") -> han_t

df %>% pivot_longer(cols = colnames(df),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% #head()
  mutate('value' = str_split_fixed(value,"\\*",2)[,2]) %>%
  mutate(Ethnity = "Korean") ->df_t

head(df_t)
head(han_t)  

df_t %>% left_join(han_t) %>%
  group_by(HLAgene) %>% #head()
  filter(HLAgene == "HLA-B") %>% #head()
  count(value) %>% arrange(-n) %>%
  ggplot(aes(x= reorder(value,n),y=n,fill=value)) + 
  geom_bar(stat = "identity") +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank()) +
  coord_flip()  ##pDPB
  
han_t %>% filter(HLAgene =="HLA-B") %>%
  count(value) %>% filter(length(value) > 6) #% arrange(n)

df_t %>% full_join(han_t) %>% #count(Ethnity)#head() c
  group_by(Ethnity,HLAgene) %>% #head()
  count(value) %>% arrange(-n) %>% #head()
  group_by(Ethnity,HLAgene) %>%
  mutate(percentage = n/sum(n)) %>% #head()
  filter()

df_t %>% #full_join(han_t) %>% #count(Ethnity)#head() c
  group_by(Ethnity,HLAgene) %>% #head()
  count(value) %>% arrange(-n) %>% #head()
  group_by(Ethnity,HLAgene) %>%
  mutate(percentage = n/sum(n)) %>% #head()
  top_n(5, percentage) %>%
  mutate(ID = paste0(HLAgene,"*",value))-> df_t_t5
  #filter(!slice_max(percentage,n=5)) %>%  head()
df_t_t5 %>% head()
  
df_t %>% #full_join(han_t) %>% #count(Ethnity)#head() c
  group_by(Ethnity,HLAgene) %>% #head()
  count(value) %>% arrange(-n) %>% #head()
  group_by(Ethnity,HLAgene) %>%
  mutate(percentage = n/sum(n)) %>% 
  mutate(ID = paste0(HLAgene,"*",value)) %>%
  filter(!(ID %in% df_t_t5$ID)) %>%
  summarise(percentage=sum(percentage), n = sum(n)) %>%
  mutate(value = "Etc") -> df_t_t6
  
df_t_t5 %>% full_join(df_t_t6) %>% #head
  ggplot(aes(x=HLAgene,y=n,fill=value)) +
  geom_bar(stat = "identity",position = "fill") +
  theme(legend.position = "none")# +
#  facet_wrap(~HLAgene,ncol = 2)

######### HLA type frequency compare with park 2016

setwd("~/Desktop/KCDC/HLA_seq/")
ref <- readxl::read_xlsx("HLAtype.freq_park2016.xlsx",sheet = 3)
df <- readxl::read_xlsx("HLAtype.freq.xlsx")
head(ref)
head(df)
df %>% mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>% count(HLAgene)
ref %>% group_by(`HLA-Gene`) %>% count(`HLA-Gene`) #head()




df %>% mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>%
  na.omit() %>% #head()
  group_by(HLAgene) %>% #head()
  mutate(Frequency = prop.table(n)) %>% #writexl::write_xlsx("HLAtype.freq.xlsx")
  mutate(type = "NIH") %>% select(-n) %>% filter(HLAgene %in% ref2$HLAgene)-> a
a
library(ggplot2)
library(ggpmisc)

#library(ggpubr)

colnames(ref) <- c("HLAgene","value","Frequency")
ref %>%  head()
a %>% head()
table(ref$HLAgene)
ref %>% filter(HLAgene != "HLA-DQB1") %>%
  mutate(type = "Park",Frequency = Frequency/100) %>%
  filter(value %in% a$value )-> ref1
head(ref1)
a %>% filter(value %in% ref$value) %>% full_join(ref1) %>% #head()
  pivot_wider(names_from = type,values_from = Frequency) %>% #head()
  ggplot(aes(x=NIH,y=Park,color=HLAgene)) + 
  geom_point() + 
  stat_smooth(method = 'lm', se=F, color='gray') + 
  ylab(element_text("Park (2016)"))
  #geom_text(x=0.05, y=0.2, label="y=0.000626 + 0.993x",color ="black") + 
  #geom_text(x=0.05, y=0.18, label="R2 = 0.95",color ="black")
  #geom_text(y=0.2,)
#  stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#  stat_regline_equation(label.y = 0.18, aes(label = ..rr.label..))

#  geom_smooth()
#  facet_wrap(~HLAgene,nrow = 2)
  #y=0.000626 + 0.993x, R2 = 0.95
library(ggVennDiagram)
head(a)
ref %>% filter(HLAgene != "HLA-DQB1") %>%
  mutate(type = "Park",Frequency = Frequency/100) -> ref2
head(ref2)

b <- a %>% full_join(ref2) 
ggVennDiagram(list(A= a$value,B=ref2$value)) + 
  theme(legend.position = "none")
head(b)

library(VennDiagram)

install.packages("ggvenn")
library(ggvenn)
head(a)
table(a$HLAgene)
head(ref2)
table(ref2$HLAgene)
ggvenn(
  list(NIH = a$value,Park =ref2$value), 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)


####
han %>% mutate(HLAtype_A.1 = paste0("A*",str_sub(HLAtype_A.1,1,-3),":",str_sub(HLAtype_A.1,-2,-1))) %>%
  mutate(HLAtype_B.1 = paste0("B*",str_sub(HLAtype_B.1,1,-3),":",str_sub(HLAtype_B.1,-2,-1))) %>%
  mutate(HLAtype_C.1 = paste0("C*",str_sub(HLAtype_C.1,1,-3),":",str_sub(HLAtype_C.1,-2,-1))) %>%
  mutate(HLAtype_DRB1.1 = paste0("DRB1*",str_sub(HLAtype_DRB1.1,1,-3),":",str_sub(HLAtype_DRB1.1,-2,-1))) %>%
  mutate(HLAtype_DPA1.1 = paste0("DPA1*",str_sub(HLAtype_DPA1.1,1,-3),":",str_sub(HLAtype_DPA1.1,-2,-1))) %>%
  mutate(HLAtype_DPB1.1 = paste0("DPB1*",str_sub(HLAtype_DPB1.1,1,-3),":",str_sub(HLAtype_DPB1.1,-2,-1))) %>%
  mutate(HLAtype_DQA1.1 = paste0("DQA1*",str_sub(HLAtype_DQA1.1,1,-3),":",str_sub(HLAtype_DQA1.1,-2,-1))) %>%
  mutate(HLAtype_DQB1.1 = paste0("DQB1*",str_sub(HLAtype_DQB1.1,1,-3),":",str_sub(HLAtype_DQB1.1,-2,-1))) %>% 
  pivot_longer(cols = colnames(han),names_to = "HLAgene") %>% na.omit() -> han1


han %>% mutate(HLAtype_A.1 = paste0("A*",str_sub(HLAtype_A.1,1,-3),":",str_sub(HLAtype_A.1,-2,-1))) %>%
  mutate(HLAtype_B.1 = paste0("B*",str_sub(HLAtype_B.1,1,-3),":",str_sub(HLAtype_B.1,-2,-1))) %>%
  mutate(HLAtype_C.1 = paste0("C*",str_sub(HLAtype_C.1,1,-3),":",str_sub(HLAtype_C.1,-2,-1))) %>%
  mutate(HLAtype_DRB1.1 = paste0("DRB1*",str_sub(HLAtype_DRB1.1,1,-3),":",str_sub(HLAtype_DRB1.1,-2,-1))) %>%
  mutate(HLAtype_DPA1.1 = paste0("DPA1*",str_sub(HLAtype_DPA1.1,1,-3),":",str_sub(HLAtype_DPA1.1,-2,-1))) %>%
  mutate(HLAtype_DPB1.1 = paste0("DPB1*",str_sub(HLAtype_DPB1.1,1,-3),":",str_sub(HLAtype_DPB1.1,-2,-1))) %>%
  mutate(HLAtype_DQA1.1 = paste0("DQA1*",str_sub(HLAtype_DQA1.1,1,-3),":",str_sub(HLAtype_DQA1.1,-2,-1))) %>%
  mutate(HLAtype_DQB1.1 = paste0("DQB1*",str_sub(HLAtype_DQB1.1,1,-3),":",str_sub(HLAtype_DQB1.1,-2,-1))) %>% 
  pivot_longer(cols = colnames(han),names_to = "HLAgene") %>% na.omit() %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "HLAtype_", replacement = "HLA-")) %>% #head()
  mutate("HLAgene" = gsub(x = HLAgene, pattern = "\\.1", replacement = "")) %>% #head()
  #mutate('value' = str_split_fixed(value,"\\*",2)[,2]) %>% head()
  group_by(HLAgene) %>% #head()
  count(value) %>% #head()
  mutate(Frequency = prop.table(n)) %>% #head()
  arrange(HLAgene,-n) %>% #head()
  mutate(type = "Han Chinese") -> hanfreq
  
df %>% mutate("HLAgene" = gsub(x = HLAgene, pattern = "NGS_", replacement = "HLA-")) %>%
  na.omit() %>% #head()
  group_by(HLAgene) %>% #head()
  mutate(Frequency = prop.table(n)) %>%
  mutate(type = "NIH") %>% select(-n) -> a



ggvenn(
  list(Korean = df$value,`Han Chinese` =han1$value), 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
table(df$HLAgene)
head(a)
table(han1$HLAgene)
head(han_t)
head(han)
head(hanfreq)
head(a)
table(a$HLAgene)

ggvenn(
  list(Korean = a$value,`Han Chinese` =hanfreq$value), 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)


hanfreq %>% select(-n) %>% filter(value %in% a$value) %>% full_join(a %>% filter(value %in% hanfreq$value) %>% mutate(type = "Korean")) %>% 
  pivot_wider(names_from = type,values_from = Frequency) %>% #head()
  ggplot(aes(x=Korean,y=`Han Chinese`,color=HLAgene)) + 
  #ggplot(aes(x=NIH,y=`Han Chinese`)) + 
  geom_point() + 
  stat_smooth(method = 'lm', se=F, color='gray') + 
  xlim(0,0.5) + 
  #xlab(element_blank())
  #ylab(element_text(" (2016)"))
  geom_text(x=0.05, y=0.4, label="y=0.0017 + 0.95x",color ="black") +
  geom_text(x=0.05, y=0.38, label="R2 = 0.84",color ="black")
#geom_text(y=0.2,)
#  stat_regline_equation(label.y = 0.2, aes(label = ..eq.label..)) +
#  stat_regline_equation(label.y = 0.18, aes(label = ..rr.label..))

#  geom_smooth()
#  facet_wrap(~HLAgene,nrow = 2)
#y=0.0017 + 0.95x, R2 = 0.84

############################################

setwd("~/Desktop/KCDC/HLA_seq/")
#writexl::write_xlsx(hanfreq,"HLAtype.freq_HAN.xlsx")
ref <- readxl::read_xlsx("HLAtype.freq_park2016.xlsx",sheet = 3)
library(corrr)

df <- readxl::read_xlsx("HLAtype.freq.xlsx")
han <- readxl::read_xlsx("HLAtype.freq_HAN.xlsx")
head(df)
head(ref)
head(han)
head(df)

df %>% filter(HLAgene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1")) %>%
  mutate(type = "NIH") %>% #head()
  filter((value %in% ref$value)) -> df1

ref %>% filter(HLAgene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1")) %>%
  mutate(type = "Park") %>% mutate(Frequency = Frequency/100) %>%
  filter((value %in% df$value)) -> ref1

df1 %>% select(-n)%>% rbind(ref1) %>% group_by(type) %>% #summary() #head()
  pivot_wider(names_from = type,values_from = Frequency) %>%
  correlate()

df %>% filter(HLAgene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1")) %>%
  mutate(type = "NIH") %>% #head()
  filter(!(value %in% ref$value)) -> df1
df1 %>% filter(Frequency > 0.05)
ref %>% filter(HLAgene %in% c("HLA-A","HLA-B","HLA-C","HLA-DRB1")) %>%
  mutate(type = "Park") %>% mutate(Frequency = Frequency/100) %>%
  filter(!(value %in% df$value)) -> ref1
head(ref1)
head(df1)
df1 %>% filter(Frequency > 0.03)
ref1 %>% filter(Frequency > 0.03)

df1 %>% select(-n)%>% rbind(ref1) %>% group_by(type) %>% #summary() #head()
  ggplot(aes(x=Frequency,fill=HLAgene)) + 
  geom_histogram(position = "dodge")  +
  facet_grid(~type)


head(han)
head(df)

head(df1)
table(han$HLAgene)
table(df$HLAgene)
df %>% mutate(type = "NIH") %>% #head()
  filter(!(value %in% han$value)) -> df1

han %>% filter(!(value %in% df$value)) -> han1


df1 %>% rbind(han1) %>% group_by(type) %>% #summary() #head()
  ggplot(aes(x=Frequency,fill=HLAgene)) + 
  geom_histogram(position = "dodge")  +
  facet_grid(~type)


df %>% mutate(type = "NIH") %>% #head()
  filter(!(value %in% han$value)) %>% arrange(-Frequency) %>% head()

han %>% filter(!(value %in% df$value))%>% 
  arrange(-Frequency) %>% head()


df %>% mutate(type = "NIH") %>% #head()
  filter((value %in% han$value)) -> df1
han %>% filter((value %in% df$value)) -> han1
head(df1)
head(han1)
df1 %>% rbind(han1) %>%select(-n) %>% group_by(type) %>% #summary() #head()
  pivot_wider(names_from = type,values_from = Frequency) %>%
  correlate()


