###############R STUDIO 단축키#####################
# 1. 코드실행
# ctrl + enter
# 2.  소스 저장
# ctrl + s
# 3. command 창 지우기
# Ctrl+L
# 4.  주석 처리/해제
# 해당 라인에 커서를 두고 ctrl + shift + c
# 5. 함수 또는 R 소스파일의 내용보기
# 확인하려는 함수 또는 R 소스파일에 커서를 올리고 F2
# 6. 실행중인 명령어 중지
# ESC
# 7. 이전 실행 명령어 창에 띄우기
# UP / DOWN
# 8. 이전 실행 명령어 확인
# CTRL + UP

#################### 패키지 설치##########################

# R 3.5 버전에서 설치
source("http://bioconductor.org/biocLite.R")
biocLite("golubEsets")
install.packages("e1071")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter")


# R 3.6 버전에서 설치

if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")
BiocManager::install("genefilter")
BiocManager::install("golubEsets")
install.packages("e1071")

#########################################################

#1. Preprocessing, feature selection
#1) Training Set & Test Set 자료 불러오기
library(golubEsets)
library(Biobase)
library(genefilter)
data(Golub_Train)
data(Golub_Test)
data(Golub_Merge)

#파일 확인하기
ls()
dim(pData(Golub_Train))
dim(pData(Golub_Test))
dim(pData(Golub_Merge))


# 1) Expression value이 100 이하이면 100로 조정하고, 16,000보다 크면 16,000으로 조정
GolubTrans <- function(eSet) {
  X <- exprs(eSet)
  X[X < 100] <- 100
  X[X > 16000] <- 16000
  X <- log2(X)
}


gTrn <- GolubTrans(Golub_Train)
gTest <- GolubTrans(Golub_Test)
gMerge <- GolubTrans(Golub_Merge)

dim(gMerge)
head(gMerge)


# 2) filtering function by Dudoit et al., 2002
# filtering by exclusion of genes with max/min<=5 or max-min<=500
mmfilt <- function(r = 5, d = 500, na.rm = TRUE) {
  function(x) {
    minval <- min(2^x, na.rm = na.rm)
    maxval <- max(2^x, na.rm = na.rm)
    (maxval/minval > r) && (maxval - minval > d)
  }
}

mmfun <- mmfilt()
ffun <- filterfun(mmfun)
sub <- genefilter(gTrn, ffun)
sum(sub)

gTrnS <- gTrn[sub, ]
gTestS <- gTest[sub, ]
gMergeS <- gMerge[sub, ]

#output 변수 y정의

Ytr <- Golub_Train$ALL.AML
Ytest <- Golub_Test$ALL.AML
Ymerge <- Golub_Merge$ALL.AML


#3) Feature selection 
#유전자마다 t-test를 하여 유의확률 0.00001보다 작은 유전자를 선택

af <- Anova(c(rep(1,27), rep(2,11)), .00001)
anovasub <- sub[sub==TRUE]
for(i in 1:sum(anovasub))
  anovasub[i] <- af(gTrnS[i,])

sum(anovasub)

gTrA <- gTrnS[anovasub, ]
gTeA <- gTestS[anovasub, ]
gMeA <- gMergeS[anovasub, ]


#4) Median Absolute Deviation (MAD) 값이 0인 유전자를 제거
whBad1 <- (apply(gTrA, 1, mad) == 0)
whBad2 <- (apply(gTeA, 1, mad) == 0)
whBad <- whBad1 | whBad2
sum(whBad) #11

gTrA <- gTrA[!whBad, ]
gTeA <- gTeA[!whBad, ]
gMeA <- gMeA[!whBad, ]

#5) Normalization
# median값을 빼서 MAD로 나누기
star <- function(x) (x - median(x))/mad(x)
TrExprs <- t(apply(gTrA, 1, star))
TeExprs <- t(apply(gTeA, 1, star))
MeExprs <- t(apply(gMeA, 1, star))



#2. Classification
#1) KNN
library(class)
knn1 <- knn(t(TrExprs), t(TeExprs), pData(Golub_Train)$ALL.AML, k = 1)
table(knn1, Golub_Test$ALL.AML)
knn3 <- knn(t(TrExprs), t(TeExprs), pData(Golub_Train)$ALL.AML, k = 3)
table(knn3, Golub_Test$ALL.AML)
knn5 <- knn(t(TrExprs), t(TeExprs), pData(Golub_Train)$ALL.AML, k = 5)
table(knn5, Golub_Test$ALL.AML)

knn3.cvpreds <- knn.cv(t(MeExprs), Ymerge,k=3)

table(knn3.cvpreds, Ymerge)

knn5.cvpreds <- knn.cv(t(MeExprs), pData(Golub_Merge)$ALL.AML,k=5)
table(knn5.cvpreds, Ymerge)


#2)svm
library(e1071)
model <- svm(t(gTrA), Golub_Train$ALL.AML, type = "C-classification", kernel="linear")
trpred <- predict(model, t(gTrA))
sum(trpred != Golub_Train$ALL.AML)
table(trpred, Golub_Train$ALL.AML)

#training dataset을 10 fold cross validation 하여 나은 classifier를 찾rl
cv_model <- svm(t(gTrA), Golub_Train$ALL.AML, type = "C-classification", kernel="linear", cross=10)
summary(cv_model)

# test dataset 분류 확
tepred <- predict(cv_model, t(gTeA))
sum(tepred != Golub_Test$ALL.AML)
table(tepred, Golub_Test$ALL.AML)


#Linear Discriminant Analysis (LDA)
library(MASS)
gTr.lda <- lda(t(TrExprs), Ytr)

# for plotting 에러 해결
par("mar")
par(mar=c(1,1,1,1))


plot(gTr.lda)

preds.lda <- predict(gTr.lda, t(TeExprs))
table(preds.lda$class, Ytest)

