### 고려대학원 2021-2 : ㄷ데이터시각화 과제
library(tidyverse)
library(readxl)

setwd("~/Desktop/KU/2021_Fall/biovisual/hw/ScienceDirect_files_19Oct2021_06-39-57.618/")


df <- read_excel("1-s2.0-S0092867421008576-mmc1.xlsx",sheet = 2)
head(df)
colnames(df)
df <- as.data.frame(df)



df %>% select(Type,Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  group_by(Type,Stage) %>%
  summarise()
    









df %>% filter(Type == "Tumor") %>% # stage 간 NAT는 NA
  select(Ethnicity,Gender,Smoking.History_modified,Stage) %>%
  #count(Ethnicity) 
  #count(Gender)
  #count(Smoking.History_modified)
  #count(Stage)
  mutate_at('Ethnicity',str_replace_all,c(han = 'asian',tssdidnotcollectthisinformation = 'NA')) %>%
  mutate_at('Ethnicity',str_replace,"white\\(caucasian\\)","caucasian") %>%
  #count(Ethnicity)
  mutate_at("Stage",str_replace_all,c(IA = 'I',IB = 'I',IIA = 'II',IIB = 'II',IIIA = 'III',IIIB = 'III')) %>%
  #count(Stage)
  gather() %>%
  filter(value != "NA") %>% 
  group_by(key,value) %>%
  summarise(n = n()) %>%
  mutate(Frequency = n/sum(n)) %>% #head()
  arrange(key,-Frequency) %>% 
  #arrange(c(1:4,6,5,9,7,8,11,10,12,13)) %>%
  arrange(c(1:4,5,6,7,8,9,11,10,12,13)) %>%
  as.data.frame() %>%
  mutate(Color = c("darkred","deepskyblue3","aquamarine4","black",
                   "darkred","deepskyblue3",
                   "darkred","deepskyblue3","aquamarine4",
                   "darkred","deepskyblue3","aquamarine4","darkslateblue")) %>%
  ggplot(aes(x=key,y=Frequency,fill=Color)) + 
  geom_bar(position="fill",stat= 'identity')
  

