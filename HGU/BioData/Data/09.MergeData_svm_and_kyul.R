kyul_tr <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/result/penalized_regression_tr_sum.csv")
kyul_ts <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/result/penalized_regression_ts_sum.csv")
svm <- read.csv("D:/biodatalab/2018-1/TCGA_with_GEO/result/svm_result.csv")
colnames(kyul_tr)
colnames(kyul_ts)
col_train <- colnames(svm)
col_train <- col_train[1:6]
col_test <- colnames(svm)
col_test <- col_test[7:9]
col_idx <- col_train[1:3]

colnames(kyul_tr) <- col_train
colnames(kyul_ts) <- col_train
colnames(kyul_ts)[4:6] <- col_test

result <- merge.data.frame(kyul_tr,kyul_ts)
result <- rbind(result,svm)
result["Foundataion","Gene_selection"]
  
?cbind

write.csv(result,"D:/biodatalab/2018-1/TCGA_with_GEO/result/result(svm_and_kyul).csv",row.names = F)
