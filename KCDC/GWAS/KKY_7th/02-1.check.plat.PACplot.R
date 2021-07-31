setwd("~/Desktop/KCDC/KKY/00.sampleInfo/")

library(stringr)
cel <- read.table("cel_files.txt",header = T)
pca <- read.table("PCA.txt",header = T)
# 043140 : DANlink
# 5507... : Teragen

DNAlink <- read.table("DANlink.2020.cel.list.txt")
tera <- read.table("2020.7th.tera.cel.list.txt")
head(cel)
tail(cel)
head(pca)
head(DNAlink)
head(tera)

DNAlink$company <- "DNAlink"
tera$company <- "Teragen"

df <- rbind(DNAlink,tera)
df$plate <- str_split_fixed(df$V1,"_",6)[,4]
df$ID <- str_replace_all(str_split_fixed(df$V1,"_",6)[,6],".CEL","")
head(str_split_fixed(df$V1,"_",6)[,6])
head(str_replace_all(str_split_fixed(df$V1,"_",6)[,6],".CEL",""))

head(df)
tail(df)

dim(table(df$plate))
head(pca)
df <- merge(df,pca,by.x = 'ID',by.y ='FID')


#pdf("1stQC/plotDATA/JG.1st.PCA.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,0,0.3)
     #     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
     #,xlim=c(-0.2, 0.5), ylim=c(-0.2,1)
     , xlab="PC1", ylab="PC2", main="KKY 7th 1st PCA", cex=1.5, pch=20)
abline(v=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.1, col=rgb(1,0,0,0.5), lty=3, lwd=2)

abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.06, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.06, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(df[df$company == "Teragen",]$PC1,
       df[df$company == "Teragen",]$PC2,
       col = rgb(1,0,0,0.3), cex = 1 , pch = 20)
points(df[df$company == "DNAlink",]$PC1,
       df[df$company == "DNAlink",]$PC2,
       col = rgb(0,0,1,0.3), cex = 1 , pch = 20)

color <- c(
  rgb(1,0,0,1),
  rgb(0,0,1,1))
list <- c("Teragen","DNAlink")
#legend(x = -0.25 ,y = 0.15,list,col = color,cex = 0.7,pch = 16)
legend(x = 0.3 ,y = 0.1,list,col = color,cex = 1,pch = 16)
dev.off()




pdf("plate_check/leftside_plate_check.txt.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,0,0.3)
     #     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
     #,xlim=c(-0.2, 0.5), ylim=c(-0.2,1)
     , xlab="PC1", ylab="PC2", main="KKY 7th 1st.PCA", cex=1.5, pch=20)
abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
#abline(v=0.06, col=rgb(1,0,0,0.5), lty=3, lwd=2)
#abline(h=0.06, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
points(df[df$PC1 < -0.05,]$PC1,
       df[df$PC1 < -0.05,]$PC2,
       col = rgb(0,0,1,1), cex = 1 , pch = 20)

dev.off()

rmlist <- df[-0.05 > df$PC1 ,]
table(rmlist$plate)
as.data.frame(t(table(rmlist$plate)))

plate_check <- as.data.frame(t(table(rmlist$plate)))
ori_plate_check <- as.data.frame(t(table(df$plate)))
head(ori_plate_check)
plate_check$ori_freq <- 0

for (i in plate_check$Var2) {
  print(i)
  plate_check[plate_check$Var2 == i,"ori_freq"] <- ori_plate_check[ori_plate_check$Var2 == i,]$Freq
}
ori_plate_check[ori_plate_check$Var2 == i,]$Freq
head(plate_check)
plate_check[plate_check$Var2 == i,"Freq"]
plate_check <- plate_check[,c("Var2","Freq","ori_freq")]
colnames(plate_check) <- c("plate","outliner_freq","number.of.sample")
write.table(plate_check,"plate_check/leftside_plate_check.txt",col.names = T,row.names = F,quote = F,sep = "\t")


pdf("plate_check/leftside_plate_check.txt.pdf", height = 10, width = 10)
plot(df$PC1, df$PC2, col=rgb(0,0,0,0.3)
     #     ,xlim=c(-0.25, 0.25), ylim=c(-0.15,0.15)
     #,xlim=c(-0.2, 0.5), ylim=c(-0.2,1)
     , xlab="PC1", ylab="PC2", main="KKY 7th 1st.PCA", cex=1.5, pch=20)
abline(v=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(v=0.06, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=0.06, col=rgb(1,0,0,0.5), lty=3, lwd=2)
abline(h=-0.05, col=rgb(1,0,0,0.5), lty=3, lwd=2)

points(df[df$PC1 < -0.05 | df$PC1 > 0.06| df$PC2 < -0.05| df$PC2 > 0.06,]$PC1,
       df[df$PC1 < -0.05 | df$PC1 > 0.06| df$PC2 < -0.05| df$PC2 > 0.06,]$PC2,
       col = rgb(1,0,0,1), cex = 1 , pch = 20)

dev.off()

rmlist <- df[0.06<df$PC1 | -0.05 > df$PC1 | df$PC2 < -0.05  | 0.06 < df$PC2,]
table(rmlist$type)
