### concordance Test for long vs Short vs merged vs KBA
library(ggplot2)
library(ggbreak)
library(tidyverse) 
library(ggpubr)
library(dplyr)
library(stringr)
library(cowplot)



setwd("~/Desktop/KCDC/HLA_seq/04.concordance/")

file_list <- list.files(pattern = ".txt")

df <- NULL
for (i in file_list) {
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
head(df)

df %>% group_by(title,type) %>% na.omit() %>%
  summarise(value = mean(Val)) %>% 
  ggplot(aes(x=title,y=value,color=title)) + 
  geom_point() + 
  facet_wrap(~type,ncol = 1) +
  theme(axis.text.x = element_blank())


df %>% group_by(title,type) %>% count(title,type) %>%
  mutate("count(SNP)" = n) %>%
  select(-n) -> nsnp 


df %>% group_by(title,type) %>% na.omit() %>%
  summarise(value = mean(Val)) %>% #head()
  left_join(nsnp) %>%
  pivot_wider(names_from = type,values_from = value) ->a

df %>% group_by(title,type) %>% na.omit() %>%
  summarise(value = mean(Val)) %>% #head()
  left_join(nsnp) %>%
  pivot_wider(names_from = type,values_from = value) %>% #head()
  ggtexttable()

df %>% group_by(title,type) %>% #na.omit() %>%
  summarise(value = mean(Val)) %>% #head()
  pivot_wider(names_from = type,values_from = value) %>%
  ggtexttable()


#writexl::write_xlsx(a,"~/Desktop/KCDC/HLA_seq/concordace.result.xlsx")

df <- readxl::read_xlsx("~/Desktop/KCDC/HLA_seq/concordace.result.xlsx",sheet = 2)
head(df)  
  

df %>% filter(Merged =="X") %>% pivot_longer(cols = c("Accuracy","Precision","Sensitivity"),values_to = "Value") %>% #head()
  mutate("VS" = str_replace_all(str_split_fixed(title,"_",2)[,2],"v"," VS ")) %>%
  mutate("title" = paste0(str_split_fixed(title,"_",2)[,1]," (",`count(SNP)`,")")) %>%
  #ggplot(aes(x=title,y=Value,color=title)) +
  ggplot(aes(x=title,y=Value,color=title)) +
  geom_point(aes(shape=VS),size = 3) + 
  facet_grid(~name) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank()) -> a

df %>% filter(Merged =="X") %>% #pivot_longer(cols = c("Accuracy","Precision","Sensitivity"),values_to = "Value") %>% #head()
  mutate("VS" = str_replace_all(str_split_fixed(title,"_",2)[,2],"v"," VS ")) %>% #head()
  #mutate("title" = paste0(str_split_fixed(title,"_",2)[,1]," (",`count(SNP)`,")")) %>% head()
  #select(1,2,10,7,8,9) %>%
  select(2,3,6,10,7,8,9) %>%
  ggtexttable() -> b

ggarrange(a,b,ncol = 1,nrow = 2,heights = c(4,4))


df %>% filter(Merged =="O") %>% pivot_longer(cols = c("Accuracy","Precision","Sensitivity"),values_to = "Value") %>% #head()
  mutate("VS" = str_replace_all(str_split_fixed(title,"_",2)[,2],"v"," VS ")) %>%
  mutate("title" = paste0(str_split_fixed(title,"_",2)[,1]," (",`count(SNP)`,")")) %>%
  #ggplot(aes(x=title,y=Value,color=title)) +
  ggplot(aes(x=title,y=Value,color=title)) +
  geom_point(aes(shape=VS),size = 3) + 
  facet_grid(~name) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank()) -> a

df %>% filter(Merged =="O") %>% #pivot_longer(cols = c("Accuracy","Precision","Sensitivity"),values_to = "Value") %>% #head()
  mutate("VS" = str_replace_all(str_split_fixed(title,"_",2)[,2],"v"," VS ")) %>% #head()
  #mutate("title" = paste0(str_split_fixed(title,"_",2)[,1]," (",`count(SNP)`,")")) %>% #head()
  select(2,3,6,10,7,8,9) %>%
  ggtexttable() -> b
b

ggarrange(a,b,ncol = 1,nrow = 2,heights = c(2,3))
