setwd("c:/Users/user/Desktop/KCDC/transplantation/re2ndQC/")
kr <- read.table("missing/KR.test.missing",header = T)
lr <- read.table("missing/LR.test.missing.missing",header = T)
head(kr)

boxplot(kr$P,lr$P)




library()

#F_missing_a : case
#F_missing_a : control

df <- subset(kr,kr$P < 1e-5)
df1 <- subset(df,df$F_MISS_A > 0.01 & df$F_MISS_U > 0.01)
df <- subset(kr,kr$P < 1e-5 & (kr$F_MISS_A > 0.01 & kr$F_MISS_U > 0.01))

head(df)


boxplot(df$F_MISS_A,df$F_MISS_U)
boxplot(kr$F_MISS_A,kr$F_MISS_U)


df <- subset(lr,lr$P < 1e-5)
df1 <- subset(df,df$F_MISS_A > 0.01 & df$F_MISS_U > 0.01)
df2 <- subset(lr,lr$P < 1e-5 & (lr$F_MISS_A > 0.01 & lr$F_MISS_U > 0.01))
