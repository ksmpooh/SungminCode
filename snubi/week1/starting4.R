setwd("C:/Users/user/git/SungminCode/snubi/")
library(MASS)

chisq.test(c(24,16),p = c(0.7,0.3))

countTable <- matrix(c(10845,189,10933,104),nrow = 2, byrow = TRUE)
countTable

rownames(countTable) <- c("Placebo","Aspirin")
colnames(countTable)<- c("No heart Attack","Heart Attack")

dim(countTable)

chisq.test(countTable)
chisq.test(countTable)$expected  # 기대 빈도를 확인하는 것


### 독립성 검정 예시2
str(birthwt)

birthwt$smoke <- factor(birthwt$smoke,label = c("Non Smoker","Smoker"))
birthwt$low <-  factor(birthwt$low,label = c("No","Yes") )

smoke_low_tb <- table(birthwt$smoke,birthwt$low)
smoke_low_tb
chisq.test(smoke_low_tb)


## fisher's exact test

TeaTasing <- matrix(c(3,1,1,3), nrow = 2)
colnames(TeaTasing) <- c("true_Milk","true_Tea")
rownames(TeaTasing) <-c("pred_Milk","pred_Tea")

chisq.test(TeaTasing)
chisq.test(TeaTasing)$expected

fisher.test(TeaTasing)


## cochran-armitage trend test

#CATT

prop.trend.test(c(13,7,21),c(42,14,28))
c(13,7,21)/c(42,14,28)
prop.trend.test(c(13,10,10),c(42,14,28))
c(13,10,10)/c(42,14,28)
chisq.test(matrix(c(13,29,10,4,10,18),ncol = 3))

## mcnemar's Test

AD = matrix(c(5,15,5,7),ncol = 2)
colnames(AD) <- c("A_AfterAD","B_AfterAD")
rownames(AD) <- c("A","B")
mcnemar.test(AD)


## ANOVA  3개 이상의 변수들의 연속형에 대해서 평균 비교하는 것 
# one-way Anova 
attach(anorexia)
str(anorexia)
change <- Postwt - Prewt
boxplot(change~Treat,col = rainbow(3))
aov.out <- aov(change~Treat)
summary(aov.out)

# 사후 검정 (post -hoc analysis)
TukeyHSD(aov.out)
plot(TukeyHSD(aov.out))

getwd()
list.files()
teaching_time <- read.table("teaching_time.txt",header = TRUE,sep = " ")
#aov.out <- aov(days~ageGroup +method,)

#반복이 잇는 이원분산 분석
str(ToothGrowth)
ToothGrowth$dose <- factor(ToothGrowth$dose)
str(ToothGrowth)
head(ToothGrowth)
aov.out <- aov(len ~supp +dose + supp:dose, data = ToothGrowth)
summary(aov.out)
