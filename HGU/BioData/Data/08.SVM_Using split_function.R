#SVM Using split function
# 1/2 train, 1/2 test
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
CV_test <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_CV_4000.csv",header = T, sep = ",")
Mean_test <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_Mean_4000.csv",header = T, sep = ",")
Var_test <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_VAR_4000.csv",header = T, sep = ",")
Annotated_308_test<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_foundation_308.csv.csv",header = T,sep = ",")
Annotated_2267_test<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_foundation_2267.csv",header = T, sep = ',')

CV_train <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_CV_4000.csv",header = T, sep = ",")
Mean_train <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_Mean_4000.csv",header = T, sep = ",")
Var_train <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_VAR_4000.csv",header = T, sep = ",")
Annotated_308_train<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_foundation_308.csv.csv",header = T,sep = ",")
Annotated_2267_train<-read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_train_foundation_2267.csv",header = T, sep = ',')
##local
#CV <- read.csv("/home/tjahn/Data/TCGA_with_GEO/train/TCGA_with_GEO_input_ensemble_test_CV_4000.csv",header = T, sep = ",")
#Mean <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_Mean_4000.csv",header = T, sep = ",")
#Var <- read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_VAR_4000.csv",header = T, sep = ",")
Annotated_308_train<-read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/TCGA_with_GEO_input_ensemble_train_foundation_308.csv",header = T,sep = ",")
#Annotated_2267<-read.csv("D:/biodatalab/2018-1/GEO_ensemble_input_final/GEO_input_ensemble_foundation_2267.csv",header = T, sep = ',')
Annotated_308_test<-read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/TCGA_with_GEO_input_ensemble_test_foundation_308.csv",header = T,sep = ",")


lists <- c(500,1000,1500,2000,2500,3000,3500,4000)
for(list in lists){
  
}
colnames(Annotated_308)
Annotated_308$index

result <- data.frame()
data <- Annotated_308
data$result <-as.factor(data$result)


for(i in 0:4){
  result <- data.frame()
  data <- read.csv(paste0("/home/tjahn/Data/input_ensemble/selected_model_",i,"_data.csv"), sep = ",", header = T)
  data$result <-as.factor(data$result)
  for(j in 1:5){
    test <-data[data$index == j,]
    train <-data[data$index != j, ]
    test <- subset(test, select = -c(index,patient,cancer_code))  
    train <- subset(train, select = -c(index,patient,cancer_code))
    #kernel 4types ????
    svm_model <- svm(result~.,data = train, kernel = "radial", cost = 1,coef.0 = 0.1 ,epsilon = 0.1)
    pred <- predict(svm_model,test)
    result_table<- table(pred,test$result)
    auc <- sum(result_table[1,1],result_table[2,2])/sum(result_table)
    #confusionMatrix(train$result,predict(svm_model))
    #svm_tune<-tune.svm(result~.,data = train,kernel = 'sigmoid',gamma = c(0.1,0.5,1,1.5,2,3,5),coef0 = c(0.1,0.5,1,2,3,5),cost = c(0.001,0.01,0.1,0.5,1,2,5,10))
    
    result[j,1] = auc
    result[j,2] = j
  }
  write.csv(result,paste0("/home/tjahn/tf_save_data/sungmin/result/SVM/radial_result_",i,".csv"))
}
