setwd("~/Desktop/KCDC/KKY/01.1stQC.PCA/PCAplot/DNAlink/")
pca <- read.table("DNAlink.PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="DNA_type0_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()



par(mfrow=c(2,4))

setwd("~/Desktop/KCDC/KKY/01.1stQC.PCA/PCAplot/all/")
pca <- read.table("PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="ALL_type0_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("t1_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="ALL_type1_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("t2_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="ALL_type2_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("t3_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="ALL_type3_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)



setwd("~/Desktop/KCDC/KKY/01.1stQC.PCA/PCAplot/merge/")
pca <- read.table("mergePCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="merge_type0_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("t1_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="merge_type1_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("t2_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="merge_type2_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("t3_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="merge_type3_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)
