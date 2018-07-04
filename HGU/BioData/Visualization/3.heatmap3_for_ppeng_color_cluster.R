library(heatmap3)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(dendextend)
library(plotrix)

## Start From Here ##
fresult <-function(result){
  if(result == 1) {"#CC0000"}
  else{"#00FF00"}
}

df<-read.csv("D:/biodatalab/2018-1/heatmap_ppeong/heatmap_ppeong_model3.csv",row.names = 1)
df <-df[order(-df$original_prob),]
df$prob_group<-cut(df$original_prob,breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),labels = c("0~10%","10~20%","20~30%","30~40%","40~50%","50~60%","60~70%","70~80%","80~90%","90~100%"),include.lowest = TRUE)
df$prob_group<-as.factor(df$prob_group)

model<-df

gc <- factor(model$cancer_code)
cancer_code<-model$cancer_code
CancerCode_color <- rainbow(19)[as.integer(gc)]
colfunc <- colorRampPalette(c("white", "darkblue"))
prob_group_col<-colfunc(10)[as.integer(model$prob_group)]


result_color<-unlist(lapply(model$real_y,fresult))

#myCols <- cbind(result_color,CancerCode_color)
myCols<- cbind(result_color,prob_group_col)
myCols<-cbind(myCols,CancerCode_color)
colnames(myCols)[1] <- "Result"
colnames(myCols)[2] <- "Prediction"
colnames(myCols)[3]<- "CancerCode"

model<-subset(model,select = -c(GSM_ID,cancer_code,real_y,original_prob,prob_group))
col_breaks = c(seq(-1,-0.0005,length=1001),
               seq(-0.000499999999,-0.0000000000001,length=1000),
               seq(0.00000000000001,0.000499999999,length=1000), 
               seq(0.0005,1,length=1001))
input<-data.matrix(model)

## remove GSM ID that has 0 sd ##
temp_data<-df[which(apply(input,1,function(x){sd(x)==0})),]
input<-input[-which(apply(input,1,function(x){sd(x)==0})),]

#png("genes_2500_by_diff.png",width = 1000, height = 1000)
### Cluster Analysis Modifying ####

#hcluster<-hclust(dist(t(input)))
#hcluster2<-hclust(dist(input))

#dend1 <- as.dendrogram(hcluster)
#dend2 <- as.dendrogram(hcluster2)

#dend1<-reorder(dend1,rowMeans(t(input),na.rm = TRUE))
#dend2<-reorder(dend2,rowMeans(input,na.rm = TRUE))
#cols_branches <- rainbow(7)
#dend1 <- color_branches(dend1, k = 7, col = cols_branches)
#dend2 <- color_branches(dend2, k = 7, col = cols_branches)

###################
#hclustfunc <- function(x) hclust(x, method="complete")
#distfunc <- function(x) dist(x, method="euclidean")

#hcluster<-hclustfunc(distfunc(t(input)))
#hcluster2<-hclustfunc(distfunc(input))
#mycl <- cutree(hcluster, k=7)
#clusterCols <- rainbow(length(unique(mycl)))
#myClusterSideBar <- clusterCols[mycl]

##################
#dend1 <- as.dendrogram(hcluster)
#dend2 <- as.dendrogram(hcluster2)

#dend1<-reorder(dend1,rowMeans(t(input),na.rm = TRUE))
#dend2<-reorder(dend2,rowMeans(input,na.rm = TRUE))
#cols_branches <- rainbow(7)
#dend1 <- color_branches(dend1, k = 7, col = cols_branches)
#dend2 <- color_branches(dend2, k = 7, col = cols_branches)


hMap<-heatmap3(t(input),
               col = bluered(4001),
               main = "Model3",
               #Colv = dend2,
               #Rowv = dend1,
               breaks = col_breaks,
               ColSideColors = myCols,
               #RowSideColors= myClusterSideBar,
               keep.dendro = T,
               
               
               scale = "none")

dend1 <- hMap$Rowv
dend2 <- hMap$Colv

#cols_branches <- rainbow(7)
#cols_branches_row <- rainbow(4)

dend1 <- color_branches(dend1, k = 15,groupLabels = F)
dend2 <- color_branches(dend2, k = 7,groupLabels = F)
#get_leaves_branches_attr(dend1)
#labels_colors(dend1) <- get_leaves_branches_col(dend1)
#plot(dend1)
#labels_colors(dend1)

#library(plotrix)
#sapply(labels_colors(dend1),color.id)


pdf("Model3.pdf")
hMap1<-heatmap3(t(input),
                col = bluered(4001),
                main = "Model3",
                Colv = dend2,
                Rowv = dend1,
                cexRow = 0.1,
                cexCol = 0.01,  
                breaks = col_breaks,
                ColSideColors = myCols,
                #RowSideColors= myClusterSideBar,
                keep.dendro = T,
                scale = "none")

dev.off()

sapply(get_leaves_branches_col(dend1),color.id)
get_leaves_branches_col(dend1)

table(get_leaves_branches_col(dend1))
labels_colors(dend1) <- get_leaves_branches_col(dend1)


??heatmap3
plot(dend1)
row_gene<-labels_colors(dend1)
table(labels_colors(dend1))
row_gene <- row_gene[which(row_gene == "#C33DBA")]
row_gene
#row_gene <- sapply(labels_colors(dend1),color.id)
#sapply("#C33DBA",color.id)



get_leaves_branches_col(dend2)
table(get_leaves_branches_col(dend2)) 
labels_colors(dend2) <- get_leaves_branches_col(dend2)
col_gene<-labels_colors(dend2)

col2 <- col_gene[which(col_gene == "#5D8400")]
col3 <- col_gene[which(col_gene == "#A86B00")]
col2 <-as.data.frame(col2)
col3 <-as.data.frame(col3)

df$cancer_code
col2<-df[row.names(col2),c("GSM_ID","cancer_code")]
col3<-df[row.names(col3),c("GSM_ID","cancer_code")]

col2
col3

col <- rbind(col2,col3)
col

write.csv(col,"d:/biodatalab/2018-1/Result/rusult_of_white_area_for_model3_in_heatmap3.csv")

df[df$GSM_ID=="GSM493919","cancer_code"]

