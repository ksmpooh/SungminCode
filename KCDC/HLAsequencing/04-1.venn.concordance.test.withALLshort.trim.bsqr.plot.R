############### 20221122 longDV short GATK hard filtering
library(ggplot2)
library(ggbreak)
library(tidyverse) 
library(ggpubr)
library(dplyr)



setwd("/Volumes/DATA/HLA_seq/05.concordance/concordance_result_20221129/")
setwd("~/")


file_list <- list.files(pattern = ".txt")
file_list
data <- read.table(file_list[1],header = T)
head(data)
#data <- read.table("QulityMetricx_forKOGO.txt", header = T)


data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))

library(tidyverse)
data$Sensitivity <- data$TP/(data$TP+data$FN)
data$Precision <- data$TP/(data$TP+data$FP)
data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)


data %>% #inner_join(maf1) %>% 
  select(Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=type,y=Val,fill=type)) +
  geom_boxplot()

data %>% #inner_join(maf1) %>% 
  select(Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
  mutate(title = "hi")

str_remove(str_remove(file_list[1],".txt"),"QulityMetrix_")


data %>% #inner_join(maf1) %>%  
  select(Sensitivity,Precision,Accuracy) %>% #head()
  summarise(Sensitivity=mean(Sensitivity))

#pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% head()
group_by(type) %>%
  summarise(mean=mean(Val))


df <- data %>% #inner_join(maf1) %>% 
  select(pos,Sensitivity,Precision,Accuracy) %>% 
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
  mutate(title = str_remove(str_remove(file_list[1],".txt"),"QulityMetrix_"))
head(df)
file_list[-1]

for (i in file_list[-1]) {
  data <- read.table(i,header = T)
  data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
  data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
  data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))
  
  data$Sensitivity <- data$TP/(data$TP+data$FN)
  data$Precision <- data$TP/(data$TP+data$FP)
  data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
  
  a <- data %>% #inner_join(maf1) %>% 
    select(pos,Sensitivity,Precision,Accuracy) %>% 
    pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>% 
    mutate(title = str_remove(str_remove(i,".txt"),"QulityMetrix_"))
  df <- rbind(df,a)
}
table(df$title)

df %>% na.omit() %>% filter(grepl("^KBA",title)) %>% #select(pos) %>% unique()#head()
  group_by(title,type) %>%
  summarise(mean = mean(Val)) %>% #head()
  ggplot(aes(x=type,y=mean,color=title,group=title)) +
  geom_point() + geom_line() +
  guides(color = guide_legend(ncol = 3,byrow = TRUE)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  scale_y_continuous(limits=c(0.98, 1)) + 
  #theme(legend.margin=margin()) %>%
  labs(title = "HLA seq. vs KBA concordance Test\n (# of SNP : 4,958)") + 
  xlab(element_blank()) + ylab(element_blank()) -> p1

#guides(fill=guide_legend(ncol=3)) #-> p1
  
df %>% na.omit() %>% filter(grepl("^KBA",title)) %>% #select(pos) %>% unique()#head()
  group_by(title,type) %>%
  summarise(mean = mean(Val)) %>% 
  pivot_wider(names_from = type,values_from = mean) %>%
  ggtexttable() -> p2
  
ggarrange(p1,p2,ncol = 1,nrow = 2,heights = c(8,3))
  

'
"KBA.LongDV"                                "KBA.shortDV"                               "KBA.shortGATK"                            
"KBA.shortGATKtrim"                         "KBA.shortGATKtrimbqsr"                     "LongDV.shortDV_onlyKBAintersect"          
"LongDV.shortDV"                            "LongDV.shortGATK_onlyKBAintersect"         "LongDV.shortGATK"                         
"LongDV.shortGATKtrim_onlyKBAintersect"     "LongDV.shortGATKtrim"                      "LongDV.shortGATKtrimbqsr_onlyKBAintersect"
"LongDV.shortGATKtrimbqsr"   
'
target <- read.table("~/Desktop/KCDC/????????????????????????????????????/??????????????????????????????/HLAseq/HLAseq.target.allPOS_PM50.txt")
head(target)
df %>% filter(grepl("KBA",title)) %>% select(pos) %>% unique()#head()
df %>% filter(!grepl("KBA",title)) %>% select(pos) %>% unique()#head()



df %>% filter(grepl("onlyKBA",title)) %>% select(pos) %>% unique()#head()

df %>% na.omit() %>% filter(grepl("onlyKBA",title)) %>% #select(pos) %>% unique()#head()
  mutate(title = str_replace(title, "_onlyKBAintersect", "")) %>% #head()
  group_by(title,type) %>% #head()
  summarise(mean = mean(Val)) %>% #head()
  ggplot(aes(x=type,y=mean,color=title,group=title)) +
  geom_point() + geom_line() +
  guides(color = guide_legend(ncol = 2,byrow = TRUE)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  scale_y_continuous(limits=c(0.98, 1)) + 
  #theme(legend.margin=margin()) %>%
  labs(title = "Long-read vs Short-read concordance Test\n (# of SNP : 4,958)") + 
  xlab(element_blank()) + ylab(element_blank()) -> p1

#guides(fill=guide_legend(ncol=3)) #-> p1

df %>% na.omit() %>% filter(grepl("onlyKBA",title)) %>% #select(pos) %>% unique()#head()
  mutate(title = str_replace(title, "_onlyKBAintersect", "")) %>% #head()
  group_by(title,type) %>%
  summarise(mean = mean(Val)) %>% 
  pivot_wider(names_from = type,values_from = mean) %>%
  ggtexttable() -> p2

ggarrange(p1,p2,ncol = 1,nrow = 2,heights = c(8,3))



##### long vs short
target <- read.table("~/Desktop/KCDC/????????????????????????????????????/??????????????????????????????/HLAseq/HLAseq.target.allPOS_PM50.txt")
df %>% filter(!grepl("KBA",title)) %>% select(pos) %>% unique()#head()

df %>% na.omit() %>% filter(!grepl("KBA",title)) %>% #select(pos) %>% unique()#head()
  #mutate(title = str_replace(title, "_onlyKBAintersect", "")) %>% #head()
  group_by(title,type) %>% #head()
  summarise(mean = mean(Val)) %>% #head()
  ggplot(aes(x=type,y=mean,color=title,group=title)) +
  geom_point() + geom_line() +
  guides(color = guide_legend(ncol = 2,byrow = TRUE)) + 
  scale_y_continuous(limits=c(0.94, 1)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  #theme(legend.margin=margin()) %>%
  labs(title = "Long-read vs Short-read Concordance Test\n (# of SNP : 35,754)") + 
  xlab(element_blank()) + ylab(element_blank()) -> p1

#guides(fill=guide_legend(ncol=3)) #-> p1

df %>% na.omit() %>% filter(!grepl("KBA",title)) %>% #select(pos) %>% unique()#head()
  #mutate(title = str_replace(title, "_onlyKBAintersect", "")) %>% #head()
  group_by(title,type) %>%
  summarise(mean = mean(Val)) %>% 
  pivot_wider(names_from = type,values_from = mean) %>%
  ggtexttable() -> p2

a1<-ggarrange(p1,p2,ncol = 1,nrow = 2,heights = c(8,3))

df %>% filter(!grepl("KBA",title)) %>% select(pos) %>% filter(pos %in% target$V1) %>% unique()#head()
df %>% na.omit() %>% filter(!grepl("KBA",title)) %>% #select(pos) %>% unique()#head()
  filter(pos %in% target$V1) %>%
  group_by(title,type) %>% #head()
  summarise(mean = mean(Val)) %>% #head()
  ggplot(aes(x=type,y=mean,color=title,group=title)) +
  geom_point() + geom_line() +
  guides(color = guide_legend(ncol = 2,byrow = TRUE)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  scale_y_continuous(limits=c(0.94, 1)) + 
  #theme(legend.margin=margin()) %>%
  labs(title = "Long-read vs Short-read Concordance Test (on Target ¡¾ 50)\n (# of SNP : 31,758)") + 
  xlab(element_blank()) + ylab(element_blank()) -> p1

#guides(fill=guide_legend(ncol=3)) #-> p1

df %>% na.omit() %>% filter(!grepl("KBA",title)) %>% #select(pos) %>% unique()#head()
  #mutate(title = str_replace(title, "_onlyKBAintersect", "")) %>% #head()
  group_by(title,type) %>%
  summarise(mean = mean(Val)) %>% 
  pivot_wider(names_from = type,values_from = mean) %>%
  ggtexttable() -> p2

a2 <- ggarrange(p1,p2,ncol = 1,nrow = 2,heights = c(8,3))


ggarrange(a1,a2,ncol = 2,nrow = 1)

#####################################
target <- read.table("~/Desktop/KCDC/????????????????????????????????????/??????????????????????????????/HLAseq/HLAseq.target.allPOS_PM50.txt")
head(target)
df %>% filter(grepl("KBA",title)) %>% select(pos) %>% unique()#head()
df %>% filter(!grepl("KBA",title)) %>% select(pos) %>% unique()#head()
df %>% mutate(title = recode(title,"KBA.LongDV"="KBA vs Long","KBA.shortGATKtrim" ="KBA vs Short",
                             "LongDV.shortGATKtrim_onlyKBAintersect" = "Long vs Short",
                             "LongDV.shortGATKtrim" = "Long vs Short (ALL)")) %>%
  group_by(title,type) %>% #filter(pos %in% target$V1) %>% #dim()#head()
  na.omit() %>%
  summarise(mean=mean(Val)) %>% #head()
  ggplot(aes(x=type,y=mean,color=title,group=title)) +
  geom_point() + geom_line() +
  scale_y_continuous(limits=c(0.91, 1)) + 
  guides(color = guide_legend(ncol = 4,byrow = TRUE)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  #theme(legend.margin=margin()) %>%
  labs(title = "Concordance Test Result\nHLA Seq. vs KBA : 5,374\nLong-read vs Short-read : 54,485") + 
  xlab(element_blank()) + ylab(element_blank())-> p1


#guides(fill=guide_legend(ncol=3)) #-> p1
table(df$title)
df %>% mutate(title = recode(title,"KBA.LongDV"="KBA vs Long","KBA.shortGATKtrim" ="KBA vs Short",
                             "LongDV.shortGATKtrim_onlyKBAintersect" = "Long vs Short",
                             "LongDV.shortGATKtrim" = "Long vs Short (ALL)")) %>%
  group_by(title,type) %>% #filter(pos %in% target$V1) %>% #dim()#head()
  na.omit() %>%
  summarise(mean=mean(Val)) %>% #head()
  pivot_wider(names_from = type,values_from = mean) %>%
  ggtexttable() -> p2

a1 <- ggarrange(p1,p2,ncol = 1,nrow = 2,heights = c(8,3))
a1
df
df %>% filter(!grepl("KBA",title)) %>% select(pos) %>% filter(pos %in% target$V1) %>% unique()#head()
df %>% filter(grepl("KBA",title)) %>% select(pos) %>% filter(pos %in% target$V1) %>% unique()#head()

df %>% mutate(title = recode(title,"KBA.LongDV"="KBA vs Long","KBA.shortGATKtrim" ="KBA vs Short",
                             "LongDV.shortGATKtrim_onlyKBAintersect" = "Long vs Short",
                             "LongDV.shortGATKtrim" = "Long vs Short (ALL)")) %>%
  group_by(title,type) %>% filter(pos %in% target$V1) %>% #dim()#head()
  na.omit() %>%
  summarise(mean=mean(Val)) %>% #head()
  ggplot(aes(x=type,y=mean,color=title,group=title)) +
  geom_point() + geom_line() +
  scale_y_continuous(limits=c(0.91, 1)) + 
  guides(color = guide_legend(ncol = 4,byrow = TRUE)) + 
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  #theme(legend.margin=margin()) %>%
  labs(title = "Concordance Test Result (on Target  ¡¾ 50) \nHLA Seq. vs KBA : 4,870\nLong-read vs Short-read : 43,081") + 
  xlab(element_blank()) + ylab(element_blank()) -> p1


#guides(fill=guide_legend(ncol=3)) #-> p1

df %>% mutate(title = recode(title,"KBA.LongDV"="KBA vs Long","KBA.shortGATKtrim" ="KBA vs Short",
                             "LongDV.shortGATKtrim_onlyKBAintersect" = "Long vs Short",
                             "LongDV.shortGATKtrim" = "Long vs Short (ALL)")) %>%
  group_by(title,type) %>% filter(pos %in% target$V1) %>% #dim()#head()
  na.omit() %>%
  summarise(mean=mean(Val)) %>% #head()
  pivot_wider(names_from = type,values_from = mean) %>%
  ggtexttable() -> p2

a2 <- ggarrange(p1,p2,ncol = 1,nrow = 2,heights = c(8,3))
a2
a1
ggarrange(a1,a2,ncol = 2,nrow = 1,widths = c(5,5))
