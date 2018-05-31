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
CancerCode_color <- rainbow(100)[as.integer(gc)]
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
png("pre_3.png",width = 2000, height = 2000)
#heatmep3 model_3, w = 1000, height = 100
#??heatmap3
heatmap3(t(input),col = bluered(100),cexRow = 0.5, cexCol = 0.5,
         main = "ppeng_3",
         margins = c(2,16),#Colv = NA, #breaks,
         #breaks = seq(-1,1,length.out= 1001),
         #ColSideColors = myCols,
         scale = "row"
) 

heatmap(t(input),col = bluered(111111))
#heatmap3
legend("topright",legend = c("Cancer","Normal"),fill = c("red","green")
       ,border = FALSE,bty = "n", y.intersp =2.0,cex = 2.0)

dev.off()

graphics.off()
 #par("mar") 
#par(mar=c(1,1,1,1))
#intersect(colnames(model_0),colnames(model_1))
#colnames(model_1)
summary(model_3)
min(model)
a<-input[1:10,1:10]
