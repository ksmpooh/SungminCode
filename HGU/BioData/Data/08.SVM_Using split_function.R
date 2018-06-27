#SVM Using split function
# 1/2 train, 1/2 test

library(e1071)
library(caret)
#################################################### VarianceTest ####################################################

GetVar<-function(genes){
  VAR<-apply(genes,2,sd)
  return(VAR)
}

TopVar <- function(genes, feature){
  VAR <- GetVar(genes)
  gene_ch <- rbind(genes, VAR)
  gene_ch_VAR <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_VAR[-nrow(gene_ch_VAR), 1:feature]
  return(gene_sub)
}

#################################################### DiffTest ####################################################
GetDiff <- function(genes, result){
  negative <- apply(genes[result==0,],2,mean)
  positive <- apply(genes[result==1,],2,mean)
  Diff <- abs(positive - negative)
  return(Diff)
}

TopDiff <- function(genes, result, feature){
  Diff <- GetDiff(genes, result)
  gene_ch <- rbind(genes, Diff)
  gene_ch_Diff <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_Diff[-nrow(gene_ch_Diff), 1:feature]
  return(gene_sub)
}

#################################################### Coefficient of Variance ##########################################
GetCV <- function(genes){
  CV<-apply(genes, 2, sd)/apply(genes, 2, mean)
  return(CV)
}

TopCV <- function(genes, feature){
  CV <- GetCV(genes)
  gene_ch <- rbind(genes, CV)
  gene_ch_CV <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_CV[-nrow(gene_ch_CV), 1:feature]
  return(gene_sub)
}

#################################################### VarianceTest ####################################################
GetMean<-function(genes){
  MEAN<-apply(genes,2,mean)
  return(MEAN)
}

TopMean <- function(genes, feature){
  MEAN <- GetMean(genes)
  gene_ch <- rbind(genes, MEAN)
  gene_ch_MEAN <- gene_ch[,rev(order(gene_ch[nrow(gene_ch),]))]
  gene_sub <-gene_ch_MEAN[-nrow(gene_ch_MEAN), 1:feature]
  return(gene_sub)
}
#####################

##data
#CV_test <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_CV_4000.csv",header = T, sep = ",")
#Mean_test <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_Mean_4000.csv",header = T, sep = ",")
#Var_test <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_VAR_4000.csv",header = T, sep = ",")
Annotated_308_test<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_foundation_308.csv.csv",header = T,sep = ",")
Annotated_2267_test<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_foundation_2267.csv",header = T, sep = ',')

#CV_train <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_CV_4000.csv",header = T, sep = ",")
#Mean_train <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_Mean_4000.csv",header = T, sep = ",")
#Var_train <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_VAR_4000.csv",header = T, sep = ",")
Annotated_308_train<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_foundation_308.csv.csv",header = T,sep = ",")
Annotated_2267_train<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_foundation_2267.csv",header = T, sep = ',')
##local
#CV <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_CV_4000.csv",header = T, sep = ",")
#Mean <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_Mean_4000.csv",header = T, sep = ",")
#Var <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_VAR_4000.csv",header = T, sep = ",")
#Annotated_308_train<-read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/TCGA_with_GEO_input_ensemble_train_foundation_308.csv",header = T,sep = ",")
#Annotated_2267<-read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_foundation_2267.csv",header = T, sep = ',')
#Annotated_308_test<-read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/TCGA_with_GEO_input_ensemble_test_foundation_308.csv",header = T,sep = ",")


result<-data.frame()
#train_result<-data.frame()
#test_result<-data.frame()
Method <- "SVM"

test_result <- data.frame()
number_of_gene <- c(500,1000,1500,2000,2500,3000,3500,4000)
lists<-c("CV","Mean","VAR")
for(list in lists){
  for(gene in number_of_gene){

    train<- read.csv(paste0("/home/tjahn/Data/TCGA_with_GEO/",train,"/TCGA_with_GEO_input_ensemble_",train,"_",list,"_",gene,".csv"),header = T, sep = ",")
    test <- read.csv(paste0("/home/tjahn/Data/TCGA_with_GEO/",test,"/TCGA_with_GEO_input_ensemble_",test,"_",list,"_",gene,".csv"),header = T, sep = ",")
    train$result <- as.factor(train$result)
    test$result <- as.factor(test$result)
    
    svm_model <- svm(result~.,data = train, kernel = "radial", cost = 1,coef.0 = 0.1 ,epsilon = 0.1)
    
    pred_train <- predict(svm_model,train)
    pred_test <- predict(svm_model,test)
    
    result_table1 <- table(pred_train,train)
    
    Train_Sensitivity <- sensitivity(pred_train)
    Train_Specificity <- specificity(pred_train)
    Train_Accuracy <- sum(result_table1[1,1],result_table1[2,2])/sum(result_table1)
    
    Gene_selection <- list
    Gene_num <- gene
    
    result_table2 <- table(pred_train,test)
    
    Test_Sensitivity <- sensitivity(pred_test)
    Test_Specificity <- specificity(pred_test)
    Test_Accuracy <- sum(result_table2[1,1],result_table2[2,2])/sum(result_table2)
    
    df <- data.frame(Gene_selection,Gene_num,Method,Train_Accuracy,Train_Sensitivity,Train_Specificity,
                     Test_Accuracy,Test_Sensitivity,Test_Specificity)
    result <-rbind(result,df)
    
  }
}

train <- Annotated_308_train
test <- Annotated_308_test

train$result <- as.factor(train$result)
test$result <- as.factor(test$result)

svm_model <- svm(result~.,data = train, kernel = "radial", cost = 1,coef.0 = 0.1 ,epsilon = 0.1)

pred_train <- predict(svm_model,train)
pred_test <- predict(svm_model,test)

result_table1 <- table(pred_train,train)

Train_Sensitivity <- sensitivity(pred_train)
Train_Specificity <- specificity(pred_train)
Train_Accuracy <- sum(result_table1[1,1],result_table1[2,2])/sum(result_table1)

Gene_selection <- "Foundation_308"
Gene_num <- "308"

result_table2 <- table(pred_train,test)

Test_Sensitivity <- sensitivity(pred_test)
Test_Specificity <- specificity(pred_test)
Test_Accuracy <- sum(result_table2[1,1],result_table2[2,2])/sum(result_table2)

df <- data.frame(Gene_selection,Gene_num,Method,Train_Accuracy,Train_Sensitivity,Train_Specificity,
                 Test_Accuracy,Test_Sensitivity,Test_Specificity)
result <-rbind(result,df)

train <- Annotated_2267_train
test <- Annotated_2267_test

train$result <- as.factor(train$result)
test$result <- as.factor(test$result)

svm_model <- svm(result~.,data = train, kernel = "radial", cost = 1,coef.0 = 0.1 ,epsilon = 0.1)

pred_train <- predict(svm_model,train)
pred_test <- predict(svm_model,test)

result_table1 <- table(pred_train,train)

Train_Sensitivity <- sensitivity(pred_train)
Train_Specificity <- specificity(pred_train)
Train_Accuracy <- sum(result_table1[1,1],result_table1[2,2])/sum(result_table1)

Gene_selection <- "Foundation_2267"
Gene_num <- "2267"

result_table2 <- table(pred_train,test)

Test_Sensitivity <- sensitivity(pred_test)
Test_Specificity <- specificity(pred_test)
Test_Accuracy <- sum(result_table2[1,1],result_table2[2,2])/sum(result_table2)

df <- data.frame(Gene_selection,Gene_num,Method,Train_Accuracy,Train_Sensitivity,Train_Specificity,
                 Test_Accuracy,Test_Sensitivity,Test_Specificity)
result <-rbind(result,df)

write.csv(result,"/home/tjahn/tf_save_data/sungmin/result/SVM/SVM_2018_06_28/result.csv",row.names = F)
