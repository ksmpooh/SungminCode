install.packages("heatmap3")

install.packages("RColorBrewer")

install.packages("gplots")

library(heatmap3)

library(RColorBrewer)

library(gplots)

 

getwd()

setwd("/Users/sungmin/Desktop/R")

gene <-read.csv("FinalData_GSM_gene_index_result.csv",header = T,sep = ",")

CancerCode <-read.csv("GPL570_sampleinfo.txt",head = T,sep= "\t")

 

row.names(gene) <-gene$X

gene<-gene[,-1]

GSM<-row.names(gene)

###Cancer Code 

row.names(CancerCode) <-CancerCode$GSM_ID

cc<-unlist(CancerCode$CANCER_CODE)

levels(cc) <-1:length(levels(cc))

cc <-as.numeric(cc)

CancerCode<-cbind(CancerCode,cc)

write.csv(CancerCode,"CancerCode_num.csv")

########

CancerCode <-CancerCode[GSM,]

gene$CancerCode<-CancerCode$cc

gene<-gene[order(-gene$result,gene$CancerCode),]

#########heatmap random sampling data

random_sample <-gene[sample(nrow(gene),100),]

a = ncol(random_sample)

random_sample <-random_sample[,((a-100):a)]

random_sample<-random_sample[order(random_sample$result,random_sample$CancerCode),]

##############################heatmap3

a= ncol(random_sample)

fresult <-function(result){

  if(result == 1) {"#CC0000"}

  else{"#00FF00"}

}

gc <- factor(random_sample$CancerCode)

CancerCode_color <- rainbow(100)[as.integer(gc)]

result_color<-unlist(lapply(random_sample$result,fresult))

 

myCols <- cbind(result_color,CancerCode_color)

colnames(myCols)[1] <- "Result"

colnames(myCols)[2] <- "CancerCode"

 

input<-data.matrix(random_sample[,(1:98)])

x = legend("topright",legend = c("Cancer","Normal"),fill = c("red","green"))

 

 

png("heatmap.png")

heatmap3(t(input),col = bluered(10),cexRow = 0.5, cexCol = 0.5,

         margins = c(2,5),Colv = NA, breaks,col

         ColSideColors = myCols,scale = "row"

) 

legend("topright",legend = c("Cancer","Normal"),fill = c("red","green"))

dev.off()

head(random_sample)

