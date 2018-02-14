#install.packages("h2o")
#install.packages("ROCR")
#install.packages("caret")
#install.packages("e1071")
#install.packages("Hmisc")
#install.packages("Dplyr")
library(Hmisc)
library(h2o)
library(ROCR)
library(caret)
library(e1071)
library(dplyr)
data <-read.csv("Total_data.csv",header = T, sep = ",")
ob <- data
#summary(ob)
#str(ob)
#str(data)
#colnames(ob)
setwd("C:/Users/bibs-student/Desktop/last_project")
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
h2o.removeAll()
for(i in 1:5){
  SNUH_test <- SNUH[SNUH$index == i,]
  SNUH_train <- SNUH[SNUH$index != i,]
  
  Asan_test <- Asan[Asan$index == i,]
  Asan_train <- Asan[Asan$index != i,]
  
  training <- rbind(SNUH_train,Asan_train)
  testing <- rbind(SNUH_test,Asan_test)
  training <- subset(training,select = -index)
  testing <-subset(testing,select = -index)
  
  training <-training[sample(1:nrow(training)),]
  testing <-testing[sample(1:nrow(testing)),]
  #colnames(training)
  training[,1] <- h2o::as.factor(training[,1])
  trData <- h2o::as.h2o(training)
  tsData <- h2o::as.h2o(testing)
  
  number_of_obs <- ncol(training)
  train_result<- h2o.deeplearning(x = 2:number_of_obs,
                                  y = "Platinum_sensitivity_in_platinum_users1",
                                  training_frame = trData,
                                  validation_frame = tsData,
                                  hidden = c(45,40,35,30,25,20),
                                  #hidden = c(30, 30, 30),
                                  #input_dropout_ratio = 0.3,
                                  #l1 = 1e-5,
                                  l2 = 1e-5,
                                  rate = 0.005,
                                  activation = "Rectifier",
                                  max_w2 = 10,
                                  epochs = 100000,
                                  stopping_rounds = 3,
                                  stopping_metric = "misclassification",
                                  #stopping_metric = "AUC",
                                  stopping_tolerance = 0.01,
                                  export_weights_and_biases = T)
  output_test[i] <- h2o.auc(train_result, valid = T)
  output_train[i] <- h2o.auc(train_result, train = T)
  
}

h2o.weights(train_result,matrix_id = 6)

h2o.shutdown(prompt=F)
summary(output_test)
summary(output_train)
sd(output_test)
sd(output_train)

output_test
output_train
