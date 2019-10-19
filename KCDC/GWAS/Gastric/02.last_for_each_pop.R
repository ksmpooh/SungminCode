

################################################

setwd("c:/Users/user/Desktop/KCDC/Gastric/QC/plotData/")
pca <- read.table("all_merge_pca.txt",header = T)
#gnomad <- read.table("1000genome_ID.txt",header = F)
samplegnomad<- read.table("1000GP_Phase3.sample",header = T)
case<-read.table("CASE_ID.txt",header = F)
control <-read.table("CONTROL_ID.txt",header = F)

#colnames(gnomad) <- "FID"
colnames(case) <- "FID"
colnames(control) <- "FID"

case$FID <-as.factor(case$FID)
control$FID <-as.factor(control$FID)

gnomad <- subset(samplegnomad,select  = c("ID","GROUP"))
colnames(gnomad) <- c("FID","GROUP")
case$GROUP <- "CASE"
control$GROUP <- "CONTROL"

df <- rbind(gnomad,case)
df <- rbind(df,control)

df <- merge(df,pca, by = "FID")
levels(df$GROUP)
pdf("../../plot/ethnic.PCA.pdf",height =10, width = 10)
plot(df$PC1,df$PC2,col = rgb(0,0,1,0.1),xlab = "PC1",ylab = "PC2",main="ethnic_PCA",
     cex.main = 3,cex = 1,pch = 16
)
points(df[df$GROUP == "AFR",]$PC1,df[df$GROUP == "AFR",]$PC2,col = rgb(1,0,1,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "AMR",]$PC1,df[df$GROUP == "AMR",]$PC2,col = rgb(0,1,1,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "EUR",]$PC1,df[df$GROUP == "EUR",]$PC2,col = rgb(0,1,0,0.3), cex = 1 , pch = 16)

points(df[df$GROUP == "SAS",]$PC1,df[df$GROUP == "SAS",]$PC2,col = rgb(1,1,0,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "CONTROL",]$PC1,df[df$GROUP == "CONTROL",]$PC2,col = rgb(0,0,1,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "CASE",]$PC1,df[df$GROUP == "CASE",]$PC2,col = rgb(1,0,0,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "EAS",]$PC1,df[df$GROUP == "EAS",]$PC2,col = rgb(0,0,0,0.3), cex = 1 , pch = 16)

color <- c(
  rgb(0,0,1,1),
  rgb(1,0,0,1),
  rgb(1,1,0,1),
  rgb(1,0,1,1),
  rgb(0,1,1,1),
  rgb(0,0,0,1),
  rgb(0,1,0,1))
list <- c("control","case","SAS","AFR","AMR","EAS","EUR")

legend(x = 0 ,y = 0.6,list,col = color,cex = 1,pch = 16)
dev.off()


##########################################################################3

