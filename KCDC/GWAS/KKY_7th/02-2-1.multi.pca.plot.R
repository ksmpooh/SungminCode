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



par(mfrow=c(4,4))

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


## merge
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


## DNAlink
setwd("~/Desktop/KCDC/KKY/01.1stQC.PCA/PCAplot/DNAlink/")
pca <- read.table("DNAlink.PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="DNAlink_type0_KKY.1st.QC_PCA", cex=1.5, pch=16)
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
     xlab="PC1", ylab="PC2", main="DNAlink_type1_KKY.1st.QC_PCA", cex=1.5, pch=16)
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
     xlab="PC1", ylab="PC2", main="DNAlink_type2_KKY.1st.QC_PCA", cex=1.5, pch=16)
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
     xlab="PC1", ylab="PC2", main="DNAlink_type3_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

setwd("~/Desktop/KCDC/KKY/01.1stQC.PCA/PCAplot/tera/")
pca <- read.table("tera.PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="Tera_type0_KKY.1st.QC_PCA", cex=1.5, pch=16)
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
     xlab="PC1", ylab="PC2", main="Tera_type1_KKY.1st.QC_PCA", cex=1.5, pch=16)
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
     xlab="PC1", ylab="PC2", main="Tera_type2_KKY.1st.QC_PCA", cex=1.5, pch=16)
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
     xlab="PC1", ylab="PC2", main="Tera_type3_KKY.1st.QC_PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)


##### 2차 정도관리를 위한 1차 정도관리

dev.off()
par(mfrow=c(2,2))

miss <-read.table("~/Desktop/KCDC/KKY/QC/DNAlink/MISS.imiss",header = T)
het <- read.table("~/Desktop/KCDC/KKY/QC/DNAlink/HET.het", header = T)


miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

#pdf("../PDF/KKY.1st_missing-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(13,22), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="DNAlink 1st QC Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=15.5, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
#dev.off()


miss <-read.table("~/Desktop/KCDC/KKY/QC/tera/MISS.imiss",header = T)
het <- read.table("~/Desktop/KCDC/KKY/QC/tera/HET.het", header = T)


miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

#pdf("../PDF/KKY.1st_missing-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(13,22), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="Teragen 1st QC Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1.5, pch=16)
abline(v=15.5, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15.5 | 17 < lowSample$HET | 0.03 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
#dev.off()

pca <- read.table("~/Desktop/KCDC/KKY/01.1stQC.PCA/PCAplot/DNAlink/t3_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="DNAlink 1st QC type3 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)

pca <- read.table("~/Desktop/KCDC/KKY/01.1stQC.PCA/PCAplot/tera/t3_PCA.txt", header=T)
plot(pca$PC1, pca$PC2, col=rgb(0,0,1,0.3),
     #xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3),
     xlab="PC1", ylab="PC2", main="DNAlink 1st QC type3 PCA", cex=1.5, pch=16)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC1,
       pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]$PC2,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
rmList <- pca[pca$PC1 < -0.1 | 0.1 < pca$PC1 | pca$PC2 < -0.1 | 0.1 < pca$PC2,]
dim(rmList)



