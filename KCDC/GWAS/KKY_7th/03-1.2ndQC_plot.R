par(mfrow=c(1,2))
setwd("~/Desktop/KCDC/KKY/02.2ndQC/DNAlink/")

miss <-read.table("MISS.imiss",header = T)
het <- read.table("HET.het", header = T)
miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")
#pdf("../PDF/KKY.1st_missing-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(13,22), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="DNAlink 2nd QC Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=15.5, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)

setwd("~/Desktop/KCDC/KKY/02.2ndQC/Tera/")

miss <-read.table("MISS.imiss",header = T)
het <- read.table("HET.het", header = T)
miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")
#pdf("../PDF/KKY.1st_missing-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(13,22), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="Teragen 2nd QC Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=15.5, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()

rmList <- lowSample[0.03 < lowSample$F_MISS | lowSample$HET < 15.5 | 17 < lowSample$HET,]
#dim(rmList)



setwd("~/Desktop/KCDC/KKY/02.2ndQC/ALL_after1stQC/")

miss <-read.table("MISS.imiss",header = T)
het <- read.table("HET.het", header = T)
miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")
#pdf("../PDF/KKY.1st_missing-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(13,22), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="ALL 2nd QC Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=15.5, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- lowSample[0.03 < lowSample$F_MISS | lowSample$HET < 15.5 | 17 < lowSample$HET,]
write.table(rmList[,c(1:2)], "rmLQSamples.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)


#### PCA
library(stringr)
DNAlink <- read.table("~/Desktop/KCDC/KKY/00.sampleInfo/DANlink.2020.cel.list.txt")
tera <- read.table("~/Desktop/KCDC/KKY/00.sampleInfo/2020.7th.tera.cel.list.txt")
DNAlink$company <- "DNAlink"
tera$company <- "Teragen"
df <- rbind(DNAlink,tera)
df$plate <- str_split_fixed(df$V1,"_",6)[,4]
df$ID <- str_replace_all(str_split_fixed(df$V1,"_",6)[,6],".CEL","")

head(df)
tail(df)

par(mfrow=c(3,3))
setwd("~/Desktop/KCDC/KKY/02.2ndQC/DNAlink/")

pca <- read.table("PCA_t0.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     xlab="PC1", ylab="PC2", main="DNAlink 2nd QC type0 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("PCA_t3.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="DNAlink 2nd QC type3 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("PCA_t4.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="DNAlink 2nd QC type4 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

setwd("~/Desktop/KCDC/KKY/02.2ndQC/Tera/")

pca <- read.table("PCA_t0.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="Teragen 2nd QC type0 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("PCA_t3.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="Teragen 2nd QC type3 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("PCA_t4.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="Teragen 2nd QC type4 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)


setwd("~/Desktop/KCDC/KKY/02.2ndQC/Merge/")

pca <- read.table("PCA_t0.txt", header=T)
pca <- merge(df,pca,by.x = 'ID',by.y ='FID')

plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     xlab="PC1", ylab="PC2", main="Merge 2nd QC type0 PCA", cex=1.5, pch=16)

points(pca[pca$company == "DNAlink",]$PC1,
       pca[pca$company == "DNAlink",]$PC2,
       col = rgb(0,0,1,0.3), cex = 1.5 , pch = 16)

points(pca[pca$company == "Teragen",]$PC1,
       pca[pca$company == "Teragen",]$PC2,
       col = rgb(1,0,0,0.3), cex = 1.5 , pch = 16)
color <- c(
  rgb(1,0,0,1),
  rgb(0,0,1,1))
list <- c("Teragen","DNAlink")
legend("bottomright",list,col = color,cex = 1,pch = 16)
#legend(x = 0 ,y = 0.48,list,col = color,cex = 1,pch = 16)
#dev.off()
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)



pca <- read.table("PCA_t3.txt", header=T)
pca <- merge(df,pca,by.x = 'ID',by.y ='FID')
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="Merge 2nd QC type3 PCA", cex=1.5, pch=16)
points(pca[pca$company == "DNAlink",]$PC1,
       pca[pca$company == "DNAlink",]$PC2,
       col = rgb(0,0,1,0.3), cex = 1.5 , pch = 16)

points(pca[pca$company == "Teragen",]$PC1,
       pca[pca$company == "Teragen",]$PC2,
       col = rgb(1,0,0,0.3), cex = 1.5 , pch = 16)
color <- c(
  rgb(1,0,0,1),
  rgb(0,0,1,1))
list <- c("Teragen","DNAlink")
#legend("bottomright",list,col = color,cex = 1,pch = 16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)


pca <- read.table("PCA_t4.txt", header=T)
pca <- merge(df,pca,by.x = 'ID',by.y ='FID')
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="Merge 2nd QC type4 PCA", cex=1.5, pch=16)
points(pca[pca$company == "DNAlink",]$PC1,
       pca[pca$company == "DNAlink",]$PC2,
       col = rgb(0,0,1,0.3), cex = 1.5 , pch = 16)

points(pca[pca$company == "Teragen",]$PC1,
       pca[pca$company == "Teragen",]$PC2,
       col = rgb(1,0,0,0.3), cex = 1.5 , pch = 16)
color <- c(
  rgb(1,0,0,1),
  rgb(0,0,1,1))
list <- c("Teragen","DNAlink")
#legend(x = 0.2 ,y = -0.1,list,col = color,cex = 1,pch = 16)
#legend("bottomright",list,col = color,cex = 1,pch = 16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)



dev.off()



par(mfrow=c(1,3))
setwd("~/Desktop/KCDC/KKY/02.2ndQC/ALL_after1stQC/")

pca <- read.table("PCA_t0.txt", header=T)
pca <- merge(df,pca,by.x = 'ID',by.y ='FID')

plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     xlab="PC1", ylab="PC2", main="ALL 2nd QC type0 PCA", cex=1.5, pch=16)

points(pca[pca$company == "DNAlink",]$PC1,
       pca[pca$company == "DNAlink",]$PC2,
       col = rgb(0,0,1,0.3), cex = 1.5 , pch = 16)

points(pca[pca$company == "Teragen",]$PC1,
       pca[pca$company == "Teragen",]$PC2,
       col = rgb(1,0,0,0.3), cex = 1.5 , pch = 16)
color <- c(
  rgb(1,0,0,1),
  rgb(0,0,1,1))
list <- c("Teragen","DNAlink")
legend("bottomright",list,col = color,cex = 1,pch = 16)
#legend(x = 0 ,y = 0.48,list,col = color,cex = 1,pch = 16)
#dev.off()
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)



pca <- read.table("PCA_t3.txt", header=T)
pca <- merge(df,pca,by.x = 'ID',by.y ='FID')
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="ALL 2nd QC type3 PCA", cex=1.5, pch=16)
points(pca[pca$company == "DNAlink",]$PC1,
       pca[pca$company == "DNAlink",]$PC2,
       col = rgb(0,0,1,0.3), cex = 1.5 , pch = 16)

points(pca[pca$company == "Teragen",]$PC1,
       pca[pca$company == "Teragen",]$PC2,
       col = rgb(1,0,0,0.3), cex = 1.5 , pch = 16)
color <- c(
  rgb(1,0,0,1),
  rgb(0,0,1,1))
list <- c("Teragen","DNAlink")
#legend("bottomright",list,col = color,cex = 1,pch = 16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)


pca <- read.table("PCA_t4.txt", header=T)
pca <- merge(df,pca,by.x = 'ID',by.y ='FID')
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.5, 0.5), ylim=c(-0.5,0.5),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="ALL 2nd QC type4 PCA", cex=1.5, pch=16)
points(pca[pca$company == "DNAlink",]$PC1,
       pca[pca$company == "DNAlink",]$PC2,
       col = rgb(0,0,1,0.3), cex = 1.5 , pch = 16)

points(pca[pca$company == "Teragen",]$PC1,
       pca[pca$company == "Teragen",]$PC2,
       col = rgb(1,0,0,0.3), cex = 1.5 , pch = 16)
color <- c(
  rgb(1,0,0,1),
  rgb(0,0,1,1))
list <- c("Teragen","DNAlink")
#legend(x = 0.2 ,y = -0.1,list,col = color,cex = 1,pch = 16)
#legend("bottomright",list,col = color,cex = 1,pch = 16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)
