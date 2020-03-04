#���Լ� ����

library(MASS)

attach(Pima.tr)
head(Pima.tr)
str(Pima.tr)
shapiro.test(bmi) ## ���Լ� �����Ҷ� ����ϴ±⺻ �Լ� 0.05���� Ŭ��� ���Լ��� �ڴ�.

qqnorm(bmi)
qqline(bmi)

#onesample t-test
t.test(bmi,mu = 30) #��հ��� 30�̾ƴϴ� ��°��� ä���Ҽ��ִ�.

bmi.ttest <- t.test(bmi,mu = 30)
names(bmi.ttest)

bmi.ttest$p.value

t.test(bmi,mu = 30, alternative = "greater") # ��������, 30���� ũ��. �븳������ �����Ҽ� �ִ�.

#two sample t-test  ..t test�� �ַ� 2���� �ϰ� �л��� �׽�Ʈ(����) 
#��л� ����
head(Pima.tr)
var.test(bmi ~ type) # bmi�� type ���� �����ȴ����� �л���� .. '~' �з� �� �� �ִ� �ڵ�����...
t.test(bmi ~type)


# paired t-test(with Anorexia data) ������ �̷���� �� ������ ������ ���� �ϰ� ���� ��

FT <- subset(anorexia,Treat == "FT")
head(FT)
shapiro.test(FT$Prewt - FT$Postwt)
t.test(FT$Prewt, FT$Postwt, paired = TRUE)



# ������ �ڷῡ�� ������ ���� ����, �����Ͱ� ���Ժ����� ������ �ʰų� ���ü��� ������, ���� ������ ���Ҷ�
# ������ �̿��Ҷ� �ַ� ����Ѵ�..

CBT <- subset(anorexia,Treat == "CBT")
head(CBT)

shapiro.test(CBT$Prewt - CBT$Postwt) #���Լ��� ������ �ʱ� ������ ���� ����� ����� �Ѵ�.(p �� 0.05���� �۴�)

wilcox.test(CBT$Prewt,CBT$Postwt,paired = TRUE)
wilcox.test(CBT$Prewt,CBT$Postwt,paired = TRUE,exact = FALSE)


placebo <- c(7,5,6,4,12)
new_drug <- c(3,6,4,2,1)

wilcox.test(placebo,new_drug, exact = FALSE) #p value �� 0.05���� ũ�⿡ �ǹ̰� ����.

# correlation test (��� �м�)

attach(iris)
head(iris)
cor(Sepal.Length,Petal.Width)
cor.test(Sepal.Length,Petal.Width)

cor(iris[,1:4])
pairs(iris[,1:4]) ## �׸�



#����ġ ������� ������

iris.na.test <- iris[,1:4]
iris.na.test[1,1] <- NA
iris.na.test[3,2] <- NA
iris.na.test[4,3] <- NA
head(iris.na.test)

cor(iris.na.test , use = "pairwise.complete.obs")

doctorA <- c(4,1,3,2,6,5,8,7)
doctorB <- c(5,3,1,2,6,4,7,8)
cor.test(doctorA,doctorB,method = "spearman")



#ȸ�ͺм�

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