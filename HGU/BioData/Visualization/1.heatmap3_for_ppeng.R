library(heatmap3)
library(RColorBrewer)
library(gplots)
for(i in 0:4){
  #assign(paste0("model_",i),read.csv(paste0("C:/Users/sungmin/Downloads/heatmap_ppeong/model_ppeong/heatmap_ppeong_model",i,".csv")))
  df <- read.csv(paste0("C:/Users/sungmin/Downloads/heatmap_ppeong/model_ppeong/heatmap_ppeong_model",i,".csv"))
  df <- df[order(-df$real_y,df$cancer_code),]
  assign(paste0("model_",i),df)
}

fresult <-function(result){
  if(result == 1) {"#CC0000"}
  else{"#00FF00"}
}

model<-model_3
gc <- factor(model$cancer_code)
cancer_code<-model$cancer_code
CancerCode_color <- rainbow(19)[as.integer(gc)]
result_color<-unlist(lapply(model$real_y,fresult))
myCols <- cbind(result_color,CancerCode_color)
colnames(myCols)[1] <- "Result"
colnames(myCols)[2] <- "CancerCode"
model<-subset(model,select = -c(X,GSM_ID,cancer_code,real_y))
#model<-subset(model,select = -c(X,GSM_ID))
#model<-model[1:10,1:10]
#model<-round(model,digits = 4)
#model<-data.matrix(subset(model,select = -X))
#input<-data.matrix(subset(model,select = -GSM_ID))
input<-data.matrix(model)
#myclust=function(c) {hclust(c,method="average")}
#input <- round(input, digits = 3)
#x = legend("topright",legend = c("Cancer","Normal"),fill = c("red","green"))
png("pre_3.png",width = 1000, height = 1000)
#heatmep3 model_3, w = 1000, height = 100
#??heatmap3
heatmap3(t(input),col = bluered(100000),cexRow = 0.5, cexCol = 0.5,
         main = "ppeng_3",
         margins = c(3,16),Colv = NA, #breaks,
         breaks = seq(-1,1,length.out= 100001),
         ColSideColors = myCols,
         keysize = 1.5,
         scale = "none"
) 

#heatmap(t(input),col = bluered(111111))
??heatmap3
legend(title = "Result","topright",legend = c("Cancer","Normal"),fill = c("red","green")
       ,border = FALSE,bty = "n", y.intersp =1.5,cex = 1.5)
legend(title = "Cancer_code","right",legend = levels(cancer_code),fill = rainbow(19))
dev.off()

graphics.off()
 #par("mar") 
#par(mar=c(1,1,1,1))
#intersect(colnames(model_0),colnames(model_1))
#colnames(model_1)
summary(model_3)
min(model)
a<-input[1:10,1:10]
table(model_3$cancer_code)
