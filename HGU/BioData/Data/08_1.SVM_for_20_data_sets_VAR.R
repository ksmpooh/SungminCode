#SVM
##test train 20

library(e1071)
library(rpart)
library(caret)

#CV <- read.csv("/home/tjahn/Data/TCGA_with_GEO/Mother/for_Model_Generation/TCGA_with_GEO_input_ensemble_CV_3000.csv",header = T,sep = ',')
#VAR <- read.csv("/home/tjahn/Data/TCGA_with_GEO/Mother/for_Model_Generation/TCGA_with_GEO_input_ensemble_VAR_3000.csv",header = T,sep = ',')
#Mean <- read.csv("/home/tjahn/Data/TCGA_with_GEO/Mother/for_Model_Generation/TCGA_with_GEO_input_ensemble_Mean_3000.csv",header = T,sep = ',')
#foundation_308 <- read.csv("/home/tjahn/Data/TCGA_with_GEO/Mother/for_Model_Generation/TCGA_with_GEO_input_ensemble_foundation_308.csv",header = T,sep = ',')
#foundation_2267 <- read.csv("/home/tjahn/Data/TCGA_with_GEO/Mother/for_Model_Generation/TCGA_with_GEO_input_ensemble_foundation_2267.csv",header = T,sep = ',')

result<-data.frame()

Model <- "SVM"
Gene_selection <- "VAR"

train_dir <- "/home/tjahn/Data/TCGA_with_GEO/index/train/"
test_dir <- "/home/tjahn/Data/TCGA_with_GEO/index/test/"
mother <- read.csv("/home/tjahn/Data/TCGA_with_GEO/Mother/for_Model_Generation/TCGA_with_GEO_input_ensemble_VAR_3000.csv",header = T,sep = ',')

for(i in 1:20){
  train_index <- read.csv(paste0(train_dir,"Train_",i,".csv"),header = T,sep = ',')
  test_index <- read.csv(paste0(test_dir,"Test_",i,".csv"),header = T,sep = ',')
  colnames(train_index) <- c("a","b")
  colnames(test_index) <- c("a","b")
  
  train <- mother[train_index$b,]
  test <- mother[test_index$b,]
  
  train$result <- as.factor(train$result)
  test$result <- as.factor(test$result)
  
  test <- subset(test, select = -c(index,patient,cancer_code))  
  train <- subset(train, select = -c(index,patient,cancer_code))
  
  svm_model <- svm(result~.,data = train, kernel = "radial", cost = 1,coef.0 = 0.1 ,epsilon = 0.1)
  
  pred_train <- predict(svm_model,train)
  pred_test <- predict(svm_model,test)
  
  result_table1 <- table(pred_train,train$result)
  
  Train_Sensitivity <- sensitivity(result_table1)
  Train_Specificity <- specificity(result_table1)
  Train_Accuracy <- sum(result_table1[1,1],result_table1[2,2])/sum(result_table1)
  
  result_table2 <- table(pred_test,test$result)
  
  Test_Sensitivity <- sensitivity(result_table2)
  Test_Specificity <- specificity(result_table2)
  Test_Accuracy <- sum(result_table2[1,1],result_table2[2,2])/sum(result_table2)
  
  index <- i
  df <- data.frame(index,Gene_selection,Model,Train_Accuracy,Train_Sensitivity,Train_Specificity,
                   Test_Accuracy,Test_Sensitivity,Test_Specificity)
  result <-rbind(result,df)
  
}


write.csv(result,"/home/tjahn/tf_save_data/sungmin/result/SVM/SVM_2018_07_05/VAR_result_20_data_set.csv",row.names = F)
