
for(i in c("0","1","2","3","4")){
  assign(paste0("model_",i),read.csv(paste0("D:/biodatalab/2018-1/input_ensemble/result_",i,".csv"),header = T,sep = ","))
}
result <-data.frame()

model_0 <-model_0[,2:3]
model_1 <-model_1[,2:3]
model_2 <-model_2[,2:3]
model_3 <-model_3[,2:3]
model_4 <-model_4[,2:3]

colnames(model_0)<-c("AUC","index")
colnames(model_1)<-c("AUC","index")
colnames(model_2)<-c("AUC","index")
colnames(model_3)<-c("AUC","index")
colnames(model_4)<-c("AUC","index")

result <-rbind(model_0,model_1)
result <-rbind(result,model_2)
result <-rbind(result,model_3)
result <-rbind(result,model_4)

result

write.csv(result,"D:/biodatalab/2018-1/input_ensemble/result.csv",row.names = T)

