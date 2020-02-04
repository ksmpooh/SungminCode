####1QC plot

setwd("c:/Users/user/Desktop/KCDC/transplantation/")

######################################



################plot
##########PCA
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
df <- df[df$type == "KR",]
pca <-read.table("1stQC/plotDATA/PCA.txt",header = T)
head(pca)

df <- merge(df,pca,by.x = "NewID",by.y ='FID')

pdf("1stQC/plot/JG.1st.PCA.onlyKR.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,1,0.3)
     #     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
    # ,xlim=c(-0.2, 0.5), ylim=c(-0.2,1)
     , xlab="PC1", ylab="PC2", main="JG.1st.onlyKR.PCA", cex=1.5, pch=20)
abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(df[df$type == "KR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "KR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(1,0,0,1), cex = 1 , pch = 20)
dev.off()



#########KR
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
df <- df[df$type == "KD",]
pca <-read.table("1stQC/plotDATA/PCA.txt",header = T)
head(pca)

df <- merge(df,pca,by.x = "NewID",by.y ='FID')

pdf("1stQC/plot/JG.1st.PCA.onlyKD.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,1,0.3)
     #     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
     # ,xlim=c(-0.2, 0.5), ylim=c(-0.2,1)
     , xlab="PC1", ylab="PC2", main="JG.1st.onlyKD.PCA", cex=1.5, pch=20)
abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(df[df$type == "KD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "KD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(1,0,0,1), cex = 1 , pch = 20)
dev.off()
#########LR
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
df <- df[df$type == "LR",]
pca <-read.table("1stQC/plotDATA/PCA.txt",header = T)
head(pca)

df <- merge(df,pca,by.x = "NewID",by.y ='FID')

pdf("1stQC/plot/JG.1st.PCA.onlyLR.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,1,0.3)
     #     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
     # ,xlim=c(-0.2, 0.5), ylim=c(-0.2,1)
     , xlab="PC1", ylab="PC2", main="JG.1st.onlyLR.PCA", cex=1.5, pch=20)
abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(df[df$type == "LR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "LR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(1,0,0,1), cex = 1 , pch = 20)
dev.off()
#########
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
df <- df[df$type == "LD",]
pca <-read.table("1stQC/plotDATA/PCA.txt",header = T)
head(pca)

df <- merge(df,pca,by.x = "NewID",by.y ='FID')

pdf("1stQC/plot/JG.1st.PCA.onlyLD.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,1,0.3)
     #     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
     # ,xlim=c(-0.2, 0.5), ylim=c(-0.2,1)
     , xlab="PC1", ylab="PC2", main="JG.1st.onlyLD.PCA", cex=1.5, pch=20)
abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(df[df$type == "LD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "LD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(1,0,0,1), cex = 1 , pch = 20)

dev.off()
