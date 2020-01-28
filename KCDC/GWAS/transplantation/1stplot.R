setwd("c:/Users/user/Desktop/KCDC/transplantation/")

######################################
setwd("c:/Users/user/Desktop/KCDC/transplantation/")

df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
cellfiles <- read.table("2nd/2nd_CEL_file_list.txt")
head(cellfiles)
colnames(cellfiles) <- c("cell_files","NewID")

out <- merge(df,cellfiles,by = "NewID")
head(out)

write.table(out,"sample_info/sample.info.with.typeNcellfiles.txt")
#######################################

df <- read.table("sample_info/sample.info.with.typeNcellfiles.txt",header = T)
head(df)

rmPCA <- read.table("2nd/rmPCA.txt")
head(rmPCA)
colnames(rmPCA) <- c("NewID","NewID1")

rmLQ <- read.table("2nd/rmLQSamples.txt")
head(rmLQ)
colnames(rmLQ) <- c("NewID","NewID1")



################plot
miss <- read.table("1stQC/plotDATA/JG.1st_SNPolisher.missing-het.imiss",header = T)
het <- read.table("1stQC/plotDATA/JG.1st_SNPolisher.missing-het.het",header = T)
head(het)
head(miss)
miss <- cbind(miss, CR=((1 - miss$F_MISS)*100))
het <- cbind(het, HET=((het$N.NM. - het$O.HOM.)/het$N.NM.)*100)

lowSample <- merge(miss, het, by="FID")
head(lowSample)

df <- read.table("sample_info/sample.info.with.type.txt",header = T)
head(df)
df <- merge(df,lowSample,by.x = "NewID",by.y = "FID")
head(df)
pdf("1stQC/plotDATA/JG.1st.miss-het.pdf", height = 7, width = 10)
plot(df$HET, df$F_MISS, 
     xlim=c(10,25), ylim=c(0,0.1), 
     xlab="heterozygosity rate",
     ylab="missing rate", main="Missing vs heterozygosity", col=rgb(0,0,0,0.3), cex=1.5, pch=20)
abline(v=15, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(v=17, col=rgb(1,0,0,1), lty=3, lwd=2)
abline(h=0.03, col=rgb(1,0,0,1), lty=3, lwd=2)

points(df[df$type == "KR" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$HET,
       df[df$type == "KR" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$F_MISS,
       col = rgb(1,0,0,1), cex = 1 , pch = 20)
points(df[df$type == "KD" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$HET,
       df[df$type == "KD" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$F_MISS,
       col = rgb(0,1,0,1), cex = 1 , pch = 20)
points(df[df$type == "LR" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$HET,
       df[df$type == "LR" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$F_MISS,
       col = rgb(1,0,1,1), cex = 1 , pch = 20)
points(df[df$type == "LD" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$HET,
       df[df$type == "LD" & (df$HET < 15 | 17 < df$HET | 0.03 < df$F_MISS),]$F_MISS,
       col = rgb(0,0,1,1), cex = 1 , pch = 20)

color <- c(
  rgb(0,0,0,0.5),
  rgb(1,0,0,1),
  rgb(0,1,0,1),
  rgb(1,0,1,1),
  rgb(0,0,1,1))
list <- c("Normal","KR","KD","LR","LD")
legend(x = 23 ,y = 0.1,list,col = color,cex = 0.7,pch = 16)
#help(legend)
dev.off()

rmList <- df[0.03 < df$F_MISS | df$HET < 15 | 17 < df$HET,]
dim(rmList)
table(rmList$type)

#write.table(rmList[,c(1:2)], "rmLQSamples.txt", col.names= FALSE, row.names=FALSE, sep="\t", quote=FALSE)

##########PCA
df <- read.table("sample_info/sample.info.with.type.txt",header = T)
pca <-read.table("1stQC/plotDATA/PCA.txt",header = T)
head(pca)
df <- merge(df,pca,by.x = "NewID",by.y ='FID')

pdf("1stQC/plotDATA/JG.1st.PCA.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,0,0.3)
#     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
     , xlab="PC1", ylab="PC2", main="JG.1st.PCA", cex=1.5, pch=20)
abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(df[df$type == "KR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "KR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(1,0,0,1), cex = 1 , pch = 20)

points(df[df$type == "KD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "KD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(0,1,0,1), cex = 1 , pch = 20)
points(df[df$type == "LR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "LR" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(0,1,1,1), cex = 1 , pch = 20)
points(df[df$type == "LD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC1,
       df[df$type == "LD" & (df$PC1 < -0.05 | 0.05 < df$PC1 | df$PC2 < -0.05 | 0.05 < df$PC2),]$PC2,
       col = rgb(0,0,1,1), cex = 1 , pch = 20)

color <- c(
  rgb(0,0,0,0.5),
  rgb(1,0,0,1),
  rgb(0,1,0,1),
  rgb(0,1,1,1),
  rgb(0,0,1,1))
list <- c("Normal","KR","KD","LR","LD")
#legend(x = -0.25 ,y = 0.15,list,col = color,cex = 0.7,pch = 16)
legend(x = 0.2 ,y = 0.9,list,col = color,cex = 0.7,pch = 16)
dev.off()
rmlist <- df[0.05<df$PC1 | -0.05 > df$PC1 | df$PC2 < -0.05  | 0.05 < df$PC2,]
table(rmlist$type)
