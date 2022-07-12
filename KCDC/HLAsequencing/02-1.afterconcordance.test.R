setwd("~/Desktop/KCDC/long_read/2022/VCF/onlySNP/DV/concordance/")
data <- read.table("QulityMetricx_longshort.txt", header = T)

data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))

data$Sensitivity <- data$TP/(data$TP+data$FN)
data$Precision <- data$TP/(data$TP+data$FP)
data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)
head(data)
data %>% select(sample,Sensitivity,Precision,Accuracy) %>%
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=type,y=Val)) +
  geom_boxplot()


library(tidyverse)

data %>% select(Sensitivity,Precision,Accuracy)

head(df)
par(mfrow=c(1,2))

data <- read.table("concordance_0002_0003.by_sample.txt",header = T)
data <- read.table("../../GATK/concordance/concordance_0002_0003.by_sample.txt",header = T)

head(data)
data <- cbind(data, TP=(data$REF.REF + data$ALT_1.ALT_1 + data$ALT_2.ALT_2))
data <- cbind(data, FN=(data$ALT_1.REF + data$ALT_2.REF + data$ALT_2.ALT_1))
data <- cbind(data, FP=(data$REF.ALT_1 + data$REF.ALT_2 + data$ALT_1.ALT_2))

data$Sensitivity <- data$TP/(data$TP+data$FN)
data$Precision <- data$TP/(data$TP+data$FP)
data$Accuracy <- data$TP/(data$TP+data$FP+data$FN)

head(data)
a <- data
b <- data

a$Tool <- "DV"
b$Tool <- "GATK"
df <- rbind(a,b)

df %>% select(Tool,sample,Sensitivity,Precision,Accuracy) %>%
  pivot_longer(cols = c(Sensitivity,Precision,Accuracy),names_to = 'type',values_to = "Val") %>%
  ggplot(aes(x=type,y=Val,fill=Tool)) +
  geom_boxplot()
  



a1 <- read.table("../AF/long.AF.txt")
a2 <- read.table("../AF/short.AF.txt")
head(a1)
head(a2)
colnames(a1)[5] <- "long"
colnames(a2)[5] <- "short"

df <- merge(a1[,c(2,5)],a2[,c(2,5)])
head(df)
df %>% select(V2,long,short) %>%
  #pivot_longer(cols = c(long,short),names_to = 'type',values_to = "AF") %>%
  ggplot(aes(x=long,y=short)) + 
  geom_tile()


if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggcorrplot")
library(ggcorrplot)
ggcorrplot(df[,c(2,3)])
cor(df$long,df$short)

head(df)
df$V2 <- as.factor(df$V2)

ggplot(df,aes(x=V2,y=V2,fill = cor(long,short))) + 
  geom_tile()

