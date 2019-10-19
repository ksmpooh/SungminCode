setwd("c:/Users/user/Desktop/KCDC/Gastric/PCA/")

pca <- read.table("PCA.txt",header = T)
gnomad <- read.table("1000genome_ID.txt",header = F)
samplegnomad<- read.table("1000GP_Phase3.sample",header = T)
case<-read.table("CASE_ID.txt",header = F)
control <-read.table("CONTROL_ID.txt",header = F)


colnames(gnomad) <- "name"
colnames(case) <- "name"
colnames(control) <- "name"

levels(samplegnomad$GROUP)
samplegnomad$color <- 0
#samplegnomad$color[samplegnomad$GROUP == 'SAS',]
samplegnomad[samplegnomad$GROUP=='SAS',]$color = rgb(1,1,0,0.3)
samplegnomad[samplegnomad$GROUP=='AFR',]$color = rgb(1,0,1,0.3)
samplegnomad[samplegnomad$GROUP=='AMR',]$color = rgb(0,1,1,0.3)
samplegnomad[samplegnomad$GROUP=='EAS',]$color = rgb(0,0,0,0.3)
samplegnomad[samplegnomad$GROUP=='EUR',]$color = rgb(0,1,0,0.3)

gnomad<- subset(samplegnomad,select = c("ID","color"))

case$name <-as.factor(case$name)
control$name <-as.factor(control$name)
#gnomad$color <- rgb(1,0,0,0.3)
colnames(gnomad) <- c("name","color")
case$color <- rgb(1,0,0,0.3)
control$color <- rgb(0,0,1,0.1)

df <- rbind(gnomad,case)
df <- rbind(df,control)
colnames(df) <- c("FID","color")
df <- merge(pca,df,by = "FID")
color <- c(rgb(1,0,0,0.3),
         rgb(0,0,1,0.3),
         rgb(1,1,0,0.3),
         rgb(1,0,1,0.3),
         rgb(0,1,1,0.3),
         rgb(0,0,0,0.3),
         rgb(0,1,0,0.3))
list <- c("case","control","SAS","AFR","AMR","EAS","EUR")
pdf("Last_PCA.pdf",height = 10,width = 10)
plot(df$PC1,df$PC2,col = df$color,xlab = "PC1",ylab = "PC2",main="Last_PCA"
     ,cex = 1,pch = 16
     )
legend(x = 0 ,y = 0.3,list,col = color,cex = 1,pch = 16)
dev.off()

################################################

setwd("c:/Users/user/Desktop/KCDC/Gastric/PCA/")
pca <- read.table("PCA.txt",header = T)
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
pdf("1kgp.gastric.PCA.pop.pdf",height =10, width = 10)
plot(df$PC1,df$PC2,col=rgb(0,0,1,0.1),xlab = "PC1",ylab = "PC2",main = "1KGP with KCHIP(GASTRIC) PCA",cex = 1, pch =16,cex.main = 3)
points(df[df$GROUP == "AFR",]$PC1,df[df$GROUP == "AFR",]$PC2,col = rgb(1,0,1,0.3), cex= 1,pch = 16)
points(df[df$GROUP == "AMR",]$PC1,df[df$GROUP == "AMR",]$PC2,col = rgb(0,1,1,0.3), cex= 1,pch = 16)
points(df[df$GROUP == "EAS",]$PC1,df[df$GROUP == "EAS",]$PC2,col = rgb(0,0,0,0.3), cex= 1,pch = 16)
points(df[df$GROUP == "EUR",]$PC1,df[df$GROUP == "EUR",]$PC2,col = rgb(0,1,0,0.3), cex= 1,pch = 16)
points(df[df$GROUP == "SAS",]$PC1,df[df$GROUP == "SAS",]$PC2,col = rgb(1,1,0,0.3), cex= 1,pch = 16)
points(df[df$GROUP == "CONTROL",]$PC1,df[df$GROUP == "CONTROL",]$PC2,col = rgb(0,0,1,0.1), cex= 1,pch = 16)
points(df[df$GROUP == "CASE",]$PC1,df[df$GROUP == "CASE",]$PC2,col = rgb(1,0,0,0.3), cex= 1,pch = 16)

color <- c(rgb(1,0,0,0.7),
           rgb(0,0,1,0.7),
           rgb(1,1,0,0.7),
           rgb(1,0,1,0.7),
           rgb(0,1,1,0.7),
           rgb(0,0,0,0.7),
           rgb(0,1,0,0.7))
list <- c("case","control","SAS","AFR","AMR","EAS","EUR")
legend(x = 0 ,y = 0.3,list,col = color,cex = 1,pch = 16)
dev.off()
##########################################################################3

