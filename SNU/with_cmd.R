library(Hmisc)
library(h2o)
library(ROCR)
library(caret)
library(e1071)
library(dplyr)
######

set.seed(123)
ob <- read.csv("FinalData.csv",header = T, sep = ",")
ob<-ob[,2:ncol(ob)]
###############################
ob_without_NAC <- subset(ob,select = -NAC)
ob<-ob_without_NAC
###############################

SNUH <-subset(ob,ob$Institution==1)
Asan <-subset(ob,ob$Institution==3)
SNUH <-subset(SNUH,select = -Institution)
Asan <-subset(Asan,select = -Institution)

##################################

output_test <-array(5)
output_train <- array(5)

count = nrow(SNUH)/5
count = as.integer(count)
remainder = nrow(SNUH)%%5

if(remainder!=0){
  Index <- c(rep(1:5,count),1:remainder)
}else{
  Index <-rep(1:5,count)}
SNUH$index <- Index
##########################

count = nrow(Asan)/5
count = as.integer(count)
remainder = nrow(Asan)%%5

if(remainder!=0){
  Index <- c(rep(1:5,count),5,4)
}else{
  Index <-rep(1:5,count)}
Asan$index <- Index

h2o.init(nthreads = -1)

#whole = 1, snu = 2, asan =3
cmd = 3
for(i in 1:5){
  SNUH_test <- SNUH[SNUH$index == i,]
  SNUH_train <- SNUH[SNUH$index != i,]
  
  Asan_test <- Asan[Asan$index == i,]
  Asan_train <- Asan[Asan$index != i,]
  if(cmd == 1){
    training <- rbind(SNUH_train,Asan_train)
    testing <- rbind(SNUH_test,Asan_test)
  }else if(cmd == 2){
    training <- SNUH_train
    testing <- SNUH_test
  }else{
    training <- Asan_train
    testing <- Asan_test
  }
  training <- subset(training,select = -index)
  testing <-subset(testing,select = -index)
  
  colnames(training)
  training[,1] <- h2o::as.factor(training[,1])
  trData <- h2o::as.h2o(training)
  tsData <- h2o::as.h2o(testing)
  
  number_of_obs <- ncol(training)
  train_result<- h2o.deeplearning(x = 2:number_of_obs,
                                  y = "Platinum_sensitivity_in_platinum_users1",
                                  training_frame = trData,
                                  hidden = c(30,30,30),
                                  input_dropout_ratio = 0.4,
                                  l1=1e-5,
                                  activation = "Rectifier",
                                  #nfold = 5,
                                  epochs = 200)
  
  ######
  predict_result_train <-h2o.predict(object = train_result,newdata = trData[,-1])
  predict_result_test <-h2o.predict(object = train_result,newdata = tsData[,-1])
  
  predict_result_train.df <-as.data.frame(predict_result_train)
  predict_result_test.df <-as.data.frame(predict_result_test)
  
  train_labels<-training[,1]
  test_labels<-testing[,1]
  
  pr_tr<-prediction(predict_result_train.df[,3],train_labels)
  pr_ts<-prediction(predict_result_test.df[,3],test_labels)
  ####
  #h2o.auc(h2o.performance(train_result))
  #h2o.performance(train_result,train = T)
  ####
  auc_train <- performance(pr_tr, measure = "auc")
  auc_train <- auc_train@y.values[[1]]
  #auc_train
  
  auc_test <- performance(pr_ts, measure = "auc")
  auc_test <- auc_test@y.values[[1]]
  #auc_test
  output_test[i] <- auc_test
  output_train[i] <- auc_train
  
}

h2o.shutdown(prompt=F)
summary(output_test)
summary(output_train)
sd(output_test)
sd(output_train)
summary(training)

