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

model<-model_2
gc <- factor(model$cancer_code)
cancer_code<-model$cancer_code
CancerCode_color <- rainbow(19)[as.integer(gc)]
result_color<-unlist(lapply(model$real_y,fresult))
myCols <- cbind(result_color,CancerCode_color)
colnames(myCols)[1] <- "Result"
colnames(myCols)[2] <- "CancerCode"
model<-subset(model,select = -c(X,GSM_ID,cancer_code,real_y))

mc <- colorRampPalette(c("red", "yellow" ,"skyblue" ,"blue"))(n = 4001)
col_breaks = c(seq(-1,-0.0005,length=1001),
               seq(-0.000499999999,-0.0000000000001,length=1000),
               #seq(-0.00000000000009,0.0000000000009,length = 1000),
               seq(0.00000000000001,0.000499999999,length=1000),            # for yellow
               seq(0.0005,1,length=1001))  

input<-data.matrix(model)
png("genes_2500_by_diff.png",width = 1000, height = 1000)
#heatmep3 model_3, w = 1000, height = 100
heatmap3(t(input),col = bluered(4001),cexRow = 0.5, cexCol = 0.5,
         main = "genes_2500_by_diff",
         margins = c(3,16),Colv = NA, #breaks,
         breaks = col_breaks,
         ColSideColors = myCols,
         #keysize = 1.5,
         scale = "none"
) 

#heatmap(t(input),col = bluered(111111))
#??heatmap3
legend(title = "Result","topright",legend = c("Cancer","Normal"),fill = c("red","green")
       ,border = FALSE,bty = "n", y.intersp =1.5,cex = 1.5)
legend(title = "Cancer_code","right",legend = levels(cancer_code),fill = rainbow(19))
#legend(title = "key-value","topleft",legend = c("-1","1"),fill = ()
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

lists <- c("gene_selected_by_mean_0.2",
          "genes_2500_var",
"genes_2500_by_diff",
"genes_foundation_265",
"genes_foundation_nci_2267")

i = 0

for(list in lists){
  df <- read.csv(paste0("C:/Users/sungmin/Downloads/heatmap_ppeong/model_ppeong/heatmap_ppeong_model",i,".csv"))
  df <- df[order(-df$real_y,df$cancer_code),]
  assign(model,df)
  gc <- factor(model$cancer_code)
  cancer_code<-model$cancer_code
  CancerCode_color <- rainbow(19)[as.integer(gc)]
  result_color<-unlist(lapply(model$real_y,fresult))
  myCols <- cbind(result_color,CancerCode_color)
  colnames(myCols)[1] <- "Result"
  colnames(myCols)[2] <- "CancerCode"
  model<-subset(model,select = -c(X,GSM_ID,cancer_code,real_y))
  
  mc <- colorRampPalette(c("red", "yellow" ,"skyblue" ,"blue"))(n = 4001)
  col_breaks = c(seq(-1,-0.0005,length=1001),
                 seq(-0.000499999999,-0.0000000000001,length=1000),
                 #seq(-0.00000000000009,0.0000000000009,length = 1000),
                 seq(0.00000000000001,0.000499999999,length=1000),            # for yellow
                 seq(0.0005,1,length=1001))  
  
  input<-data.matrix(model)
  png(paste0(list,".png"),width = 1000, height = 1000)
  #heatmep3 model_3, w = 1000, height = 100
  heatmap3(t(input),col = bluered(4001),cexRow = 0.5, cexCol = 0.5,
           main = list,
           margins = c(3,16),Colv = NA, #breaks,
           breaks = col_breaks,
           ColSideColors = myCols,
           #keysize = 1.5,
           scale = "none"
  ) 
  
  #heatmap(t(input),col = bluered(111111))
  #??heatmap3
  legend(title = "Result","topright",legend = c("Cancer","Normal"),fill = c("red","green")
         ,border = FALSE,bty = "n", y.intersp =1.5,cex = 1.5)
  legend(title = "Cancer_code","right",legend = levels(cancer_code),fill = rainbow(19))
  #legend(title = "key-value","topleft",legend = c("-1","1"),fill = ()
  dev.off()
  i = i+1
}
