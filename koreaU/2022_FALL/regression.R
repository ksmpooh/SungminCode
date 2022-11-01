
setwd("~/Desktop/KU/2022_FALL/genetics/regression/")
library(tidyverse)
library(readxl)

ref <-  read_excel("1-s2.0-S0092867420314008-mmc1.xlsx",sheet = 1)
ref
df <-  read_excel("1-s2.0-S0092867420314008-mmc1.xlsx",sheet = 2)
rna <-  read_excel("1-s2.0-S0092867420314008-mmc1.xlsx",sheet = 3)
wes <-  read_excel("1-s2.0-S0092867420314008-mmc1.xlsx",sheet = 4)

head(wes)
head(df)
colnames(df)
library(GGally)
install.packages("GGally")
table(df$Number.of.non.synonymous.Mutations)
colnames(df)
colnames(wes)
df %>% select()
wes %>% mutate(CASE=ifelse(Type=="Tumor",1,0)) %>%#head()
  select(c(4:14,1)) %>% #head()
  #pivot_longer(cols = 4:14,names_to = "var",values_to = "value") + #%>% head()
  ggpairs(columns=1:11,aes(color = Type,alpha = 0.5),
          upper = list(continuous = wrap("cor", size = 2.5)))

wes %>% mutate(CASE=ifelse(Type=="Tumor",1,0)) %>% #colnames()
  select(c(1,4,7,9,11,13)) %>% #head()
  #rename()
  #pivot_longer(cols = 4:14,names_to = "var",values_to = "value") + #%>% head()
  ggpairs(columns=2:6,aes(color = Type,alpha = 0.5),
          upper = list(continuous = wrap("cor", size = 3)))



head(rna)
rna %>% select(1,3,5,6) %>%#$head()
  ggpairs(2:4)

head(df)
colnames(df)
grep("*|*",df$TMT.Plex)
grepl("\\|",df$TMT.Plex)

sum(str_split(df[grepl("\\|",df$TMT.Plex),]$TMT.Plex,"\\|"))
library(naniar)
df %>% select(Sample.ID,TMT.Plex,NMF.Cluster.Membership.Score,Age.in.Month,Gender,
              ESTIMATE.TumorPurity,
              CIBERSORT.Absolute.Score,
              ESTIMATE.Immune.Score,
              xCell.Immune.Score,
              ESTIMATE.Stromal.Score,
              xCell.Stromal.Score,Number.of.non.synonymous.Mutations,
              Chromosome.INstability.index.CIN.,
              Stemness.Score) %>%
  mutate(Gender = ifelse(Gender == 'female','2',ifelse(Gender=="male",'1',NA))) %>%
  mutate(TMT.Plex = ifelse(grepl("\\|",TMT.Plex),NA,TMT.Plex)) %>% #head()
  replace_with_na_all(~.x %in% c("unknown")) %>% #head()
  mutate(Age.in.Month = as.numeric(Age.in.Month)) %>%
  ggpairs(2:14,upper = list(continuous = wrap("cor", size = 2.5)))


rna1 <- rna %>% select(1,3,5,6)
head(rna1)
colnames(rna1)[2:4]<-c("Total_reads","mapping_rates","No_of_transcripts")
df %>% select(Sample.ID,
              ESTIMATE.TumorPurity,
              #CIBERSORT.Absolute.Score,
              ESTIMATE.Immune.Score,
              xCell.Immune.Score,
              #ESTIMATE.Stromal.Score,
              #xCell.Stromal.Score,
              #Number.of.non.synonymous.Mutations,
              #Stemness.Score
              #Chromosome.INstability.index.CIN.,
              Number.of.non.synonymous.Mutations) %>%
  rename("No.of.non.syn.Mutation" = Number.of.non.synonymous.Mutations) %>%
  inner_join(rna1) %>%  #head()
  ggpairs(2:8,upper = list(continuous = wrap("cor", size = 6)))

table(df$Tumor.Stage)
df_lm<-df %>% select(Sample.ID, Tumor.Stage, Age.in.Month,#TMT.Plex,
        #df %>% select(Sample.ID, Tumor.Stage,
              ESTIMATE.TumorPurity,
              ESTIMATE.Immune.Score,
              xCell.Immune.Score,
              Number.of.non.synonymous.Mutations,
              Stemness.Score) %>% #count(Tumor.Stage) #%>%
  replace_with_na_all(~.x %in% c("unknown")) %>% #head()
  mutate(Age.in.Month = as.numeric(Age.in.Month)) %>%
  mutate(Tumor.Stage = ifelse(Tumor.Stage %in% c("Stage IA"),1,ifelse(Tumor.Stage %in% c("Stage IIA","Stage IIB"),2,ifelse(Tumor.Stage=="unknown",NA,3)))) %>%
  #mutate(TMT.Plex == ifelse(str_detect(TMT.Plex,"\\|"),NA,TMT.Plex)) %>%
  inner_join(rna1) #head()
  #rename("Total_reads" = )#head()
  #lm(Number.of.non.synonymous.Mutations ~ xCell.Immune.Score)
summary(df_lm)
apply(df_lm, 2, sd)
apply(df_lm, 2, mean)

a<-df_lm %>%
  summarise(across(4:11, mean))
a
head(df_lm)
colnames(df_lm)
model <- lm(ESTIMATE.TumorPurity ~ Total_reads+mapping_rates+No_of_transcripts,data=df_lm)
model <- lm(ESTIMATE.TumorPurity ~ scale(Total_reads)+scale(mapping_rates)+No_of_transcripts,data=df_lm)


model <- lm(Stemness.Score ~ No_of_transcripts,data=df_lm)
model <- lm(Stemness.Score ~ Number.of.non.synonymous.Mutations,data=df_lm)
model <- lm(Stemness.Score ~ No_of_transcripts + ESTIMATE.TumorPurity +ESTIMATE.Immune.Score+ xCell.Immune.Score + scale(Total_reads)+scale(mapping_rates)+Number.of.non.synonymous.Mutations,data=df_lm)

model <- lm(Stemness.Score ~ No_of_transcripts,data=df_lm)
model <- lm(Stemness.Score ~ No_of_transcripts + ESTIMATE.TumorPurity +ESTIMATE.Immune.Score+ xCell.Immune.Score + scale(Total_reads)+scale(mapping_rates)+Number.of.non.synonymous.Mutations,data=df_lm)

model
summary(model)


t <- read.delim("test.tab",sep = "\t")
