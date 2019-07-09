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
ob <- ob[!is.na(ob$Platinum_sensitivity_in_platinum_users1),]
#a <- data[data$Institution == 3,]
#a[,c("Platinum_sensitivity_in_platinum_users1","Bone.")]
a <- 0
b= 1
c= 1
for (i in 1:ncol(ob)) {
  if(sum(is.na(ob[,i]))>70){
    a[b] <- i  
    #obj[c] <- colnames(ob[i])
    c = c + 1
    b = b + 1
  }
}
a <- as.integer(a)
ob <- ob[,-a]
ob <- subset(ob,select =-c(Birth_date,Dx_date,
                           Start_Date_1st_regimen,End_Date_1st_regimen,
                           Last_FU_date,Op_date,CA125_initial,NLR,PLR,Brain.,
                           IntraPeritoneal_CTx,X1st_Regimen,Bone.,
                           Personal_history_of_gynecologic_cancer,Smoking,Personal_history_of_breast_cancer))
#summary(ob$br)
#str(ob)
#ob[is.na(ob$CA125_initial),"CA125_initial"] <- mean(data$CA125_initial,na.rm = T)
ob[is.na(ob$WBC),"WBC"] <- mean(data$WBC,na.rm = T)
ob[is.na(ob$Hb),"Hb"] <- mean(data$Hb,na.rm = T)
ob[is.na(ob$PLT),"PLT"] <- mean(data$PLT,na.rm = T)
ob[is.na(ob$Segmented_neutrophil),"Segmented_neutrophil"] <- mean(data$Segmented_neutrophil,na.rm = T)
ob[is.na(ob$Lymphocyte),"Lymphocyte"] <- mean(data$Lymphocyte,na.rm = T)
ob[is.na(ob$Monocyte),"Monocyte"] <- mean(data$Monocyte,na.rm = T)
#ob[is.na(ob$NLR),"NLR"] <- mean(data$NLR,na.rm = T)
ob[is.na(ob$MLR),"MLR"] <- mean(data$MLR,na.rm = T)
#ob[is.na(ob$PLR),"PLR"] <- mean(data$PLR,na.rm = T)
ob[is.na(ob$log_CA125initial),"log_CA125initial"] <- mean(data$log_CA125initial,na.rm = T)
ob[is.na(ob$count_Segmentedneutrophil),"count_Segmentedneutrophil"] <- mean(data$count_Segmentedneutrophil,na.rm = T)
ob[is.na(ob$count_Lymphocyte),"count_Lymphocyte"] <- mean(data$count_Lymphocyte,na.rm = T)
ob[is.na(ob$count_Monocyte),"count_Monocyte"] <- mean(data$count_Monocyte,na.rm = T)
ob[is.na(ob$logNLR),"logNLR"] <- mean(data$logNLR,na.rm = T)
ob[is.na(ob$logPLR),"logPLR"] <- mean(data$logPLR,na.rm = T)

#str(ob)
#summary(ob)

Max_Value <-function(Value){
  return (as.integer(names(Value)[which.max(Value)]))
}

#t <- table(ob$FIGO1234)
#names(t)[which.max(t)]

ob[is.na(ob$Dyslipidemia),"Dyslipidemia"] <- Max_Value(table(data$Dyslipidemia))
#ob[is.na(ob$Large_bowel_resection),"Large_bowel_resection"] <- Max_Value(table(data$Large_bowel_resection))
#ob[is.na(ob$Upper_abdominal_surgery),"Upper_abdominal_surgery"] <- Max_Value(table(data$Upper_abdominal_surgery))
ob[is.na(ob$Ov_involvement),"Ov_involvement"] <- Max_Value(table(data$Ov_involvement))
ob[is.na(ob$Ov_surface),"Ov_surface"] <- Max_Value(table(data$Ov_surface))
ob[is.na(ob$Origin),"Origin"] <- Max_Value(table(data$Origin))
ob[is.na(ob$Cytology_Ascites_or_Peritoneal_washing),"Cytology_Ascites_or_Peritoneal_washing"] <- Max_Value(table(data$Cytology_Ascites_or_Peritoneal_washing))
ob[is.na(ob$No_of_harvested_LNs),"No_of_harvested_LNs"] <- Max_Value(table(data$No_of_harvested_LNs))
ob[is.na(ob$Large_bowel_resection),"Large_bowel_resection"] <- Max_Value(table(data$Large_bowel_resection))
ob[is.na(ob$Upper_abdominal_surgery),"Upper_abdominal_surgery"] <- Max_Value(table(data$Upper_abdominal_surgery))
ob[is.na(ob$Ov_involvement),"Ov_involvement"] <- Max_Value(table(data$Ov_involvement))
ob[is.na(ob$Ov_surface),"Ov_surface"] <- Max_Value(table(data$Ov_surface))
ob[is.na(ob$Ov_capsule),"Ov_capsule"] <- Max_Value(table(data$Ov_capsule))
ob[is.na(ob$Tube),"Tube"] <- Max_Value(table(data$Tube))
ob[is.na(ob$Lung.),"Lung."] <- Max_Value(table(data$Lung.))
ob[is.na(ob$Liver_parenchyme),"Liver_parenchyme"] <- Max_Value(table(data$Liver_parenchyme))
ob[is.na(ob$Supraclavicular_LN),"Supraclavicular_LN"] <- Max_Value(table(data$Supraclavicular_LN))
ob[is.na(ob$Residual_tumor_site_1st_debulking),"Residual_tumor_site_1st_debulking"] <- Max_Value(table(data$Residual_tumor_site_1st_debulking))
ob[is.na(ob$Recurrence),"Recurrence"] <- Max_Value(table(data$Recurrence))
ob[is.na(ob$Spleen0),"Spleen0"] <- Max_Value(table(data$Spleen0))
ob[is.na(ob$Spleen01),"Spleen01"] <- Max_Value(table(data$Spleen01))
ob[is.na(ob$Histology_SEMC),"Histology_SEMC"] <- Max_Value(table(data$Histology_SEMC))
ob[is.na(ob$Histology_serous),"Histology_serous"] <- Max_Value(table(data$Histology_serous))
ob[is.na(ob$small_bowel_and_mesentery01vs2),"small_bowel_and_mesentery01vs2"] <- Max_Value(table(data$small_bowel_and_mesentery01vs2))
ob[is.na(ob$small_bowel_and_mesentery0vs12),"small_bowel_and_mesentery0vs12"] <- Max_Value(table(data$small_bowel_and_mesentery0vs12))
ob[is.na(ob$Other_colon_except_rectosigmoid01vs2),"Other_colon_except_rectosigmoid01vs2"] <- Max_Value(table(data$Other_colon_except_rectosigmoid01vs2))
ob[is.na(ob$Other_colon_except_rectosigmoid0vs12),"Other_colon_except_rectosigmoid0vs12"] <- Max_Value(table(data$Other_colon_except_rectosigmoid0vs12))
ob[is.na(ob$Liver_surface01vs2),"Liver_surface01vs2"] <- Max_Value(table(data$Liver_surface01vs2))
ob[is.na(ob$Liver_surface0vs12),"Liver_surface0vs12"] <- Max_Value(table(data$Liver_surface0vs12))
ob[is.na(ob$Residual_tumor_size_1st_debulking0),"Residual_tumor_size_1st_debulking0"] <- Max_Value(table(data$Residual_tumor_size_1st_debulking0))
ob[is.na(ob$Residual_tumor_size_1st_debulking01),"Residual_tumor_size_1st_debulking01"] <- Max_Value(table(data$Residual_tumor_size_1st_debulking01))

#summary(ob)
#str(ob)
#str(data)
#colnames(ob)
set.seed(123)
ob <- ob %>% select(Platinum_sensitivity_in_platinum_users1, everything())
write.csv(ob,"FinalData.csv")


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

for(i in 1:5){
  SNUH_test <- SNUH[SNUH$index == i,]
  SNUH_train <- SNUH[SNUH$index != i,]
  
  Asan_test <- Asan[Asan$index == i,]
  Asan_train <- Asan[Asan$index != i,]
  
  training <- rbind(SNUH_train,Asan_train)
  testing <- rbind(SNUH_test,Asan_test)
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


