#### plot
setwd("~/Desktop/KCDC/transplantation/QCrepliation_2020/1st/")
setwd("~/Desktop/KCDC/transplantation/QCrepliation_2020/2nd/")

miss <-read.table("MISS.imiss",header = T)
het <- read.table("HET.het", header = T)


miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

pdf("../PDF/JG.KR_rep.QC_SNPolisher_miss-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(13,23), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=16.3, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17.8, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 16.3 | 17.8 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 16.3 | 17.8 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()

rmList <- lowSample[0.03 < lowSample$F_MISS | lowSample$HET < 16.3 | 17.8 < lowSample$HET,]
#dim(rmList)

#write.table(rmList[,c(1:2)], "rmLQSamples.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)


## pca plot

setwd("~/Desktop/KCDC/transplantation/QCrepliation_2020/1st/")

pca <- read.table("PCA.txt", header=T)

head(pca)

pdf("../PDF/JG.1st.QC_PCA.pdf", height = 10, width = 10)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     xlab="PC1", ylab="PC2", main="2nd.QC_PCA", cex=1.5, pch=16)
abline(v=-0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()

rmList <- pca[pca$PC1 < -0.13 | 0.13 < pca$PC1 | pca$PC2 < -0.13 | 0.13 < pca$PC2,]
dim(rmList)

write.table(rmList[,c(1:2)], "rmPCA.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)

q()






### ethnic

setwd("~/Desktop/KCDC/transplantation/QCrepliation_2020/1st/ethnic/")

#pca <- read.table("PCA.txt",header = T)
pca <- read.table("JG.1st.ethnic.PCA.txt",header = T)
#gnomad <- read.table("../INPUTs/1000genome_ID.txt",header = F)
samplegnomad<- read.table("1000GP_Phase3.sample",header = T)
#case<-read.table("CASE_ID.txt",header = F)
case <- read.table("JG.intersect.fam",header = F)
case <- as.data.frame(case$V1)
#control <-read.table("CONTROL_ID.txt",header = F)
control <- as.data.frame(samplegnomad$ID)

colnames(case) <- "FID"
#colnames(control) <- "FID"

case$FID <- as.factor(case$FID)
#control$FID <- as.factor(control$FID)

gnomad <- subset(samplegnomad,select = c("ID","GROUP"))
colnames(gnomad) <- c("FID","GROUP")
case$GROUP <- "CASE"
#control$GROUP <- "CONTROL"

df <- rbind(gnomad,case)
df <- rbind(df,control)

df <- merge(pca,df,by = "FID")
pdf("../PDF/JG.2nd.KR.case.control.1kgp.ethnic.PCA.pdf",height = 10,width = 10)

plot(df$PC1,df$PC2,col = rgb(0,0,1,0.1),xlab = "PC1",ylab = "PC2",main="Ethnic PCA",
     cex.main = 3,cex = 1,pch = 16
)
points(df[df$GROUP == "AFR",]$PC1,df[df$GROUP == "AFR",]$PC2,col = rgb(0,0,1,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "AMR",]$PC1,df[df$GROUP == "AMR",]$PC2,col = rgb(1,1,0,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "EUR",]$PC1,df[df$GROUP == "EUR",]$PC2,col = rgb(1,0,1,0.3), cex = 1 , pch = 16)

points(df[df$GROUP == "SAS",]$PC1,df[df$GROUP == "SAS",]$PC2,col = rgb(0,1,0,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "EAS",]$PC1,df[df$GROUP == "EAS",]$PC2,col = rgb(0,1,1,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "CONTROL",]$PC1,df[df$GROUP == "CONTROL",]$PC2,col = rgb(0,0,0,0.5), cex = 1 , pch = 16)
points(df[df$GROUP == "CASE",]$PC1,df[df$GROUP == "CASE",]$PC2,col = rgb(1,0,0,0.5), cex = 1 , pch = 16)

color <- c(
        rgb(0,0,0,1),
        rgb(1,0,0,1),
        rgb(0,1,0,1),
        rgb(0,0,1,1),
        rgb(1,1,0,1),
        rgb(0,1,1,1),
        rgb(1,0,1,1))
list <- c("control","JG.1st","SAS","AFR","AMR","EAS","EUR")
legend(x = 0.4 ,y = 0.5,list,col = color,cex = 1,pch = 16)
dev.off()


points(df[df$FID == "NIH20KT1693",]$PC1,df[df$FID == "NIH20KT1693",]$PC2,col=rgb(0,0,0,1),cex = 1 , pch = 16)

rmList <- df[df$PC1 < -0.1 | 0.1 < df$PC1 | df$PC2 < -0.1 | 0.1 < df$PC2,]
rmList[rmList$GROUP == "CASE",]
rmList <- df[ -0.25 < df$PC1 | df$PC2 > -0.16 | 0.1 < df$PC2,]
#dim(rmList)

write.table(rmList[,c(1:2)], "rmPCA.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)







####2nd QC freq

setwd("~/Desktop/KCDC/transplantation/QCrepliation_2020/03.merge/")

case <- read.table("CASE_intersect_freq.frq",header=T)
control <- read.table("CONTROL_intersect_freq.frq",header=T)
data <- merge(control,case,by="SNP")
#png("control&case_frequency.png",height = 700,width=700)

plot(data$MAF.x,data$MAF.y,xlab = "Control",ylab = "Case",main = "Control & Case Frequency")
abline(a = 0.05,b = 1, col = 'red',lty = 2)
abline(a = -0.05,b = 1, col = 'red',lty = 2)
points(data[data$MAF.x-data$MAF.y >= 0.05 | data$MAF.x - data$MAF.y <= -0.05,]$MAF.x,
       data[data$MAF.x-data$MAF.y >= 0.05 | data$MAF.x - data$MAF.y <= -0.05,]$MAF.y,
       col = 'blue', cex = 1, pch = 1)

####2nd merge PCA
setwd("~/Desktop/KCDC/transplantation/QCrepliation_2020/03.merge/")

pca <- read.table("PCA.txt", header=T)

head(pca)

pdf("../PDF/JG.1st.QC_PCA.pdf", height = 10, width = 10)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     xlab="PC1", ylab="PC2", main="KR replication (case_control) PCA", cex=1.5, pch=16)
abline(v=-0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.13, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()

rmList <- pca[pca$PC1 < -0.13 | 0.13 < pca$PC1 | pca$PC2 < -0.13 | 0.13 < pca$PC2,]
dim(rmList)

write.table(rmList[,c(1:2)], "rmPCA.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)

q()



## ethnic
setwd("~/Desktop/KCDC/transplantation/QCrepliation_2020/03.merge/ethnic/")

pca <- read.table("PCA.txt",header = T)
#pca <- read.table("JG.1st.ethnic.PCA.txt",header = T)
#gnomad <- read.table("../INPUTs/1000genome_ID.txt",header = F)
samplegnomad<- read.table("../../01.1stQC/ethnic/1000GP_Phase3.sample",header = T)
#case<-read.table("CASE_ID.txt",header = F)
case <- read.table("JG.intersect.fam",header = F)
case <- read.table("case.txt")
head(case)
#case <- as.data.frame(case$V1)
#control <-read.table("CONTROL_ID.txt",header = F)
control <- as.data.frame(samplegnomad$ID)

colnames(case) <- "FID"
#colnames(control) <- "FID"

case$FID <- as.factor(case$FID)
#control$FID <- as.factor(control$FID)

gnomad <- subset(samplegnomad,select = c("ID","GROUP"))
colnames(gnomad) <- c("FID","GROUP")
case$GROUP <- "CASE"
#control$GROUP <- "CONTROL"

df <- rbind(gnomad,case)
df <- rbind(df,control)

df <- merge(pca,df,by = "FID")
#pdf("../PDF/JG.2nd.KR.case.control.1kgp.ethnic.PCA.pdf",height = 10,width = 10)

plot(df$PC1,df$PC2,col = rgb(0,0,1,0.1),xlab = "PC1",ylab = "PC2",main="Ethnic PCA",
     cex.main = 3,cex = 1,pch = 16
)
points(df[df$GROUP == "AFR",]$PC1,df[df$GROUP == "AFR",]$PC2,col = rgb(0,0,1,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "AMR",]$PC1,df[df$GROUP == "AMR",]$PC2,col = rgb(1,1,0,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "EUR",]$PC1,df[df$GROUP == "EUR",]$PC2,col = rgb(1,0,1,0.3), cex = 1 , pch = 16)

points(df[df$GROUP == "SAS",]$PC1,df[df$GROUP == "SAS",]$PC2,col = rgb(0,1,0,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "EAS",]$PC1,df[df$GROUP == "EAS",]$PC2,col = rgb(0,1,1,0.3), cex = 1 , pch = 16)
points(df[df$GROUP == "CONTROL",]$PC1,df[df$GROUP == "CONTROL",]$PC2,col = rgb(0,0,0,0.5), cex = 1 , pch = 16)
points(df[df$GROUP == "CASE",]$PC1,df[df$GROUP == "CASE",]$PC2,col = rgb(1,0,0,0.3), cex = 1 , pch = 16)

color <- c(
#        rgb(0,0,0,1),
        rgb(1,0,0,1),
        rgb(0,1,0,1),
        rgb(0,0,1,1),
        rgb(1,1,0,1),
        rgb(0,1,1,1),
        rgb(1,0,1,1))
list <- c("control","JG.1st","SAS","AFR","AMR","EAS","EUR")
list <- c("control","JG.1st","SAS","AFR","AMR","EAS","EUR")
legend(x = 0.4 ,y = 0.5,list,col = color,cex = 1,pch = 16)
dev.off()




rmList <- df[df$PC1 < -0.1 | 0.1 < df$PC1 | df$PC2 < -0.1 | 0.1 < df$PC2,]
#dim(rmList)

write.table(rmList[,c(1:2)], "rmPCA.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)
