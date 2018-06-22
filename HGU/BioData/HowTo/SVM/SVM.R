### SVM
library(e1071)
result <- data.frame()

data <- read.csv("D:/biodatalab/2018-1/input_ensemble/selected_model_3_data.csv", sep = ",", header = T)
data$result <-as.factor(data$result)
j = 1;
test <-data[data$index == j,]
train <-data[data$index != j, ]
test <- subset(test, select = -c(index,patient,cancer_code))  
train <- subset(train, select = -c(index,patient,cancer_code))
#kernel 4types ????
svm_model <- svm(result~.,data = train, kernel = "radial", cost = 1,coef.0 = 0.1 ,epsilon = 0.1)
pred <- predict(svm_model,test)
result_table<- table(pred,test$result)
result_table
auc <- sum(result_table[1,1],result_table[2,2])/sum(result_table)
#confusionMatrix(train$result,predict(svm_model))
#svm_tune<-tune.svm(result~.,data = train,kernel = 'sigmoid',gamma = c(0.1,0.5,1,1.5,2,3,5),coef0 = c(0.1,0.5,1,2,3,5),cost = c(0.001,0.01,0.1,0.5,1,2,5,10))

result[j,1] = auc
result[j,2] = j

