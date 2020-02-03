setwd("c:/Users/user/Desktop/KCDC/transplantation/")








pca <- read.table("ethical/pca.txt",header = T)
#gnomad <- read.table("1000genome_ID.txt",header = F)
samplegnomad<- read.table("ethical/1000GP_Phase3.sample",header = T)
case<-read.table("ethical/case_ID.txt",header = F)

casetype <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(case)
head(casetype)
colnames(case) <- "FID"

casetype <- casetype[,c("NewID","type")]
colnames(casetype) <- c("FID","GROUP")
case$FID <- as.factor(case$FID)

case<- merge(case,casetype,by = 'FID')
head(case)
kgp <- subset(samplegnomad,select = c("ID","GROUP"))
colnames(kgp) <- c("FID","GROUP")

df <- rbind(kgp,case)

df <- merge(pca,df,by = "FID")

pdf("ethical/JG.ethnic_PCA.pdf",height = 10,width = 10)
#png("ethical/JG.ethnic_PCA.png")
plot(df$PC1,df$PC2,col = rgb(0,0,0,0.3),xlab = "PC1",ylab = "PC2",main="Ethnic_PCA",
     cex.main = 3,cex = 0.8,pch = 20
)
points(df[df$GROUP == "AFR",]$PC1,df[df$GROUP == "AFR",]$PC2,col = rgb(0.8,0.3,0.5,0.3), cex = 0.8 , pch = 20)
points(df[df$GROUP == "AMR",]$PC1,df[df$GROUP == "AMR",]$PC2,col = rgb(1,0,1,0.3), cex = 0.8 , pch = 20)
points(df[df$GROUP == "EUR",]$PC1,df[df$GROUP == "EUR",]$PC2,col = rgb(1,1,0,0.3), cex = 0.8 , pch = 20)
points(df[df$GROUP == "SAS",]$PC1,df[df$GROUP == "SAS",]$PC2,col = rgb(0.5,0.8,0.3,0.3), cex = 0.8 , pch = 20)
points(df[df$GROUP == "EAS",]$PC1,df[df$GROUP == "EAS",]$PC2,col = rgb(0,0,0,0.3), cex = 0.8 , pch = 20)


points(df[df$GROUP == "KD",]$PC1,df[df$GROUP == "KD",]$PC2,col = rgb(0,1,0,0.1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "LR",]$PC1,df[df$GROUP == "LR",]$PC2,col = rgb(0,1,1,0.1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "LD",]$PC1,df[df$GROUP == "LD",]$PC2,col = rgb(0,0,1,0.1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "KR",]$PC1,df[df$GROUP == "KR",]$PC2,col = rgb(1,0,0,0.1), cex = 0.8 , pch = 20)
color <- c(
  rgb(0.8,0.3,0.5,1),
  rgb(1,0,1,1),
  rgb(1,1,0,1),
  rgb(0.5,0.8,0.3,1),
  rgb(0,0,0,1),
  rgb(1,0,0,1),
  rgb(0,1,0,1),
  rgb(0,1,1,1),
  rgb(0,0,1,1))
list <- c("AFR","AMR","EUR","SAS","EAS","KR","KD","LR","LD")

legend(x = -0.8 ,y = 0.4,list,col = color,cex = 1.3,pch = 16)
dev.off()



############################################tiny

pdf("ethical/JG.ethnic_PCA.zoom.pdf",height = 10,width = 10)
#png("ethical/JG.ethnic_PCA.png")
plot(df$PC1,df$PC2,col = rgb(0,0,0,0.3),xlab = "PC1",ylab = "PC2",main="Ethnic_PCA",
     xlim=c(0, 0.2), ylim=c(-0.1,0.1),
     cex.main = 3,cex = 0.8,pch = 20
)
points(df[df$GROUP == "AFR",]$PC1,df[df$GROUP == "AFR",]$PC2,col = rgb(0.8,0.3,0.5,1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "AMR",]$PC1,df[df$GROUP == "AMR",]$PC2,col = rgb(1,0,1,1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "EUR",]$PC1,df[df$GROUP == "EUR",]$PC2,col = rgb(1,1,0,1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "SAS",]$PC1,df[df$GROUP == "SAS",]$PC2,col = rgb(0.5,0.8,0.3,1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "EAS",]$PC1,df[df$GROUP == "EAS",]$PC2,col = rgb(0,0,0,1), cex = 0.8 , pch = 20)


points(df[df$GROUP == "KD",]$PC1,df[df$GROUP == "KD",]$PC2,col = rgb(0,1,0,1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "LR",]$PC1,df[df$GROUP == "LR",]$PC2,col = rgb(0,1,1,1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "LD",]$PC1,df[df$GROUP == "LD",]$PC2,col = rgb(0,0,1,1), cex = 0.8 , pch = 20)
points(df[df$GROUP == "KR",]$PC1,df[df$GROUP == "KR",]$PC2,col = rgb(1,0,0,1), cex = 0.8 , pch = 20)
color <- c(
#  rgb(0.8,0.3,0.5,1),
#  rgb(1,0,1,1),
#  rgb(1,1,0,1),
#  rgb(0.5,0.8,0.3,1),
  rgb(0,0,0,1),
  rgb(1,0,0,1),
  rgb(0,1,0,1),
  rgb(0,1,1,1),
  rgb(0,0,1,1))
#list <- c("AFR","AMR","EUR","SAS","EAS","KR","KD","LR","LD")
list <- c("EAS","KR","KD","LR","LD")
legend(x = 0 ,y = 0.1,list,col = color,cex = 1,pch = 16)
dev.off()

############################################miss-het
miss <- read.table("2nd/JG.2nd.QC_snpolisher_miss-het.imiss",header = T)
het <- read.table("2nd/JG.2nd.QC_snpolisher_miss-het.het",header = T)
head(het)
head(miss)
miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")

pdf("2nd/JG.2nd.miss-het.pdf", height = 7, width = 10)
plot(lowSample$HET, lowSample$F_MISS, xlim=c(10,25), ylim=c(0,0.1), xlab="heterozygosity rate",
     ylab="missing rate", main="2ndQC.Missing vs heterozygosity", col=rgb(0,0,1,0.3), cex=1, pch=20)
abline(v=15, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)
points(lowSample[lowSample$HET < 15.4 | 17 < lowSample$HET | 0.3 < lowSample$F_MISS,]$HET,
       lowSample[lowSample$HET < 15.4 | 17 < lowSample$HET | 0.3 < lowSample$F_MISS,]$F_MISS,
       col=rgb(1,0,0,0.3), cex=1.5, pch=16)
dev.off()

