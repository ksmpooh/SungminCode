library(e1071)
library(rpart)
library(caret)
data <- read.csv("d:/biodatalab/2018-1/input_ensemble/selected_model_0_data.csv", sep = ",", header = T)
data$result <-as.factor(data$result)

test <-data[data$index == 1,]
train <-data[data$index != 1, ]
test <- subset(test, select = -c(index,patient,cancer_code))  
train <- subset(train, select = -c(index,patient,cancer_code))

data <- subset(data, select = -c(index,patient,cancer_code))

svm_model <- svm(result~.,data = train, kernel = "sigmoid", cost = 0.005, gamma= 0.5,coef.0 = 0.1 ,epsilon = 0.1)
svm_model1 <- svm(result~.,data = train, kernel = "sigmoid", cost = 0.005, gamma= 0.5,coef.0 = 0.1 ,epsilon = 0.1)
svm_model2 <- svm(result~.,data = train, kernel = "sigmoid", cost = 0.005, gamma= 0.5,coef.0 = 0.1,epsilon = 0.1)
svm_model3 <- svm(result~.,data = train, kernel = "sigmoid", cost = 0.005, gamma= 0.5,coef.0 = 0.1 ,epsilon = 0.1)
svm_model4 <- svm(result~.,data = train, kernel = "sigmoid", cost = 0.005, gamma= 0.5,coef.0 = 0.1 ,epsilon = 0.1)
#best cost = 0.005,gamma 영향 없음,codf.도 영향 없음
pred <- predict(svm_model,test)
pred1 <- predict(svm_model1,test)
pred2 <- predict(svm_model2,test)
pred3 <- predict(svm_model3,test)
pred4 <- predict(svm_model4,test)
a<-confusionMatrix(train$result,predict(svm_model))
b<-confusionMatrix(train$result,predict(svm_model1))
c<-confusionMatrix(train$result,predict(svm_model2))
d<-confusionMatrix(train$result,predict(svm_model3))
e<-confusionMatrix(train$result,predict(svm_model4))

a
b
c
d
e

pred
pred1
pred2
pred3
pred4

list
#svm_train<-tune.svm(result~.,data = train,kernel = 'sigmoid',gamma = c(0.1,0.3,0.5,0.7,1,1.5,2,3,5),coef0 = c(0.1,0.3,0.5,0.7,1,1.5,2,3,5),cost = c(0.001,0.005,0.01,0.05,0.1,0.3,0.5,1,1.5,2,5,10))
#svm_data<-tune.svm(result~.,data = data,kernel = 'sigmoid',gamma = c(0.1,0.3,0.5,0.7,1,1.5,2,3,5),coef0 = c(0.1,0.3,0.5,0.7,1,1.5,2,3,5),cost = c(0.001,0.005,0.01,0.05,0.1,0.3,0.5,1,1.5,2,5,10)
library(twitteR)
ta <- twListToDF(a)
                   
                   write(svm_train$best.parameters,"/home/tjahn/tf_save_data/sungmin/svm_train.csv")
                   write(svm_data$best.parameters,"/home/tjahn/tf_save_data/sungmin/svm_data.csv")
                   
                   
svm_model_data <- svm(result~.,data = data, kernel = "sigmoid", type = "C-classification",cost = 0.005, gamma= 0.5,coef.0 = 0.1 ,epsilon = 0.1)
#e<-confusionMatrix(train$result,predict(svm_model))
#a<-predict(svm_model_data,data$result)
confusionMatrix(data$result,predict(svm_model_data))

pred <- predict(svm_model,test)
result_table<- table(pred,test$result)
result_table
confusionMatrix(test$result,pred)
