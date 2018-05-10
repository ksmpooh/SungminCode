
for(i in c("0","1","2","3","4")){
  assign(paste0("model_",i),read.csv(paste0("D:/biodatalab/2018-1/input_ensemble/Radial/Radial_result_",i,".csv"),header = T,sep = ","))
}
result <-data.frame()

#model_0 <-model_0[,2:3]
#model_1 <-model_1[,2:3]
#model_2 <-model_2[,2:3]
#model_3 <-model_3[,2:3]
#model_4 <-model_4[,2:3]

colnames(model_0)<-c("ensemble_model", "accuracy","method")
colnames(model_1)<-c("ensemble_model", "accuracy","method")
colnames(model_2)<-c("ensemble_model", "accuracy","method")
colnames(model_3)<-c("ensemble_model", "accuracy","method")
colnames(model_4)<-c("ensemble_model", "accuracy","method")

result <-rbind(model_0,model_1)
result <-rbind(result,model_2)
result <-rbind(result,model_3)
result <-rbind(result,model_4)

result

result[,"ensemble_model"]<-rep(0:4,each =5)
result[,"method"] <- rep("DNN",times=5)
write.csv(result,"D:/biodatalab/2018-1/input_ensemble/result_SVM_radial.csv",row.names = F)

