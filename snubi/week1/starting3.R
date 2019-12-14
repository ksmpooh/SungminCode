#정규성 검정

library(MASS)

attach(Pima.tr)
head(Pima.tr)
str(Pima.tr)
shapiro.test(bmi) ## 정규성 검증할때 사용하는기본 함수 0.05보다 클경우 정규성을 뛴다.

qqnorm(bmi)
qqline(bmi)

#onesample t-test
t.test(bmi,mu = 30) #평균값이 30이아니다 라는것을 채택할수있다.

bmi.ttest <- t.test(bmi,mu = 30)
names(bmi.ttest)

bmi.ttest$p.value

t.test(bmi,mu = 30, alternative = "greater") # 단측검증, 30보다 크다. 대립가설을 선택할수 있다.

#two sample t-test  ..t test는 주로 2개를 하고 분산을 테스트(검증) 
#등분산 검증
head(Pima.tr)
var.test(bmi ~ type) # bmi가 type 별로 구문된다음에 분산검증 .. '~' 분류 할 수 있는 자동검증...
t.test(bmi ~type)


# paired t-test(with Anorexia data) 쌍으로 이루어진 두 변수의 차리를 검정 하고 싶을 때

FT <- subset(anorexia,Treat == "FT")
head(FT)
shapiro.test(FT$Prewt - FT$Postwt)
t.test(FT$Prewt, FT$Postwt, paired = TRUE)



# 연속형 자료에서 비모수적 가설 검정, 데이터가 정규분포를 따르지 않거나 샘플수가 적을때, 순위 측도를 비교할때
# 순위를 이용할때 주로 사용한다..

CBT <- subset(anorexia,Treat == "CBT")
head(CBT)

shapiro.test(CBT$Prewt - CBT$Postwt) #정규성을 따르지 않기 떄문에 비모수 방법을 따라야 한다.(p 가 0.05보다 작다)

wilcox.test(CBT$Prewt,CBT$Postwt,paired = TRUE)
wilcox.test(CBT$Prewt,CBT$Postwt,paired = TRUE,exact = FALSE)


placebo <- c(7,5,6,4,12)
new_drug <- c(3,6,4,2,1)

wilcox.test(placebo,new_drug, exact = FALSE) #p value 가 0.05보다 크기에 의미가 없다.

# correlation test (상관 분석)

attach(iris)
head(iris)
cor(Sepal.Length,Petal.Width)
cor.test(Sepal.Length,Petal.Width)

cor(iris[,1:4])
pairs(iris[,1:4]) ## 그림



#결측치 있을경우 상관계수

iris.na.test <- iris[,1:4]
iris.na.test[1,1] <- NA
iris.na.test[3,2] <- NA
iris.na.test[4,3] <- NA
head(iris.na.test)

cor(iris.na.test , use = "pairwise.complete.obs")

doctorA <- c(4,1,3,2,6,5,8,7)
doctorB <- c(5,3,1,2,6,4,7,8)
cor.test(doctorA,doctorB,method = "spearman")



#회귀분석

newTest <-  c(50,55,60,65,70,75,80,85,90,95,100)
standardTest <- c(61,61,59,71,80,76,90,106,98,100,114)
dat <- data.frame(newTest, standardTest)

plot(standardTest ~newTest, data = dat, xlim = c(0,110), ylim = c(0,120))
dat.lm <- lm(standardTest~newTest, data= dat)
summary(dat.lm)

coef(dat.lm)
predict(dat.lm,newdata = data.frame(newTest = 80))

install.packages("lmtest")
library(lmtest)
dwtest(standardTest ~ newTest, data = dat)
