#Original signal intensity에서 transformed signal intensity로 변환할 수 있는 공식
#x축: Contrast = log2(SNP-A / SNP-B)
#y축: Size = (log2(SNP-A) + log2(SNP-B))/2

library(stringr)
library(tidyverse)
setwd("~/Desktop/KCDC/KKY/01.1stQC_clusterplot/")



df <- read_table("summary_for.cluster.txt")
df <- as.data.frame(df)
df2 <- as.data.frame(t(df))
colnames(df2) <- df$probeset_id
head(df2[2:nrow(df2),1:2])
df2 <- df2[2:nrow(df2),]
head(df2)
df2 <- read.table("summary_for.cluster_transposon.txt",header = T)
#### cluster plot

df <- read.table("clusterplot.data.x.y.txt",header = T)
markerID <- read.table("markerID.txt")
ref <- read.table("TEST2_raw_processing.AAABBB.txt",header = T)

head(df)[1:3]
df$ID <- str_sub(str_split_fixed(df$CEL,"_",6)[,6],1,-5)

#df <- df[df$ID %in% ref$FID,]
row.names(df)<-df$ID




png("plot/test.png",height = 600,width = 600)
plot(df[,'AX.106718326.contrast'],df[,'AX.106718326.size'], cex=0.5, pch=0,col=rgb(1,1,1,0.5),
     xlab = "Contrast = log2(SNP-A / SNP-B)",ylab = "Size = (log2(SNP-A) + log2(SNP-B))/2",
     )
#dev.off()
#plot(xlab = "Contrast = log2(SNP-A / SNP-B)",ylab = "Size = (log2(SNP-A) + log2(SNP-B))/2")
points(df[ref[ref[,i] == 0,]$FID,paste0(i,'.contrast')],
       df[ref[ref[,i] == 0,]$FID,paste0(i,'.size')],
          cex=0.5, pch=16,col=rgb(0,1,0,0.5))
points(df[ref[ref[,i] == 1,]$FID,paste0(i,'.contrast')],
       df[ref[ref[,i] == 1,]$FID,paste0(i,'.size')],
       cex=0.5, pch=16,col=rgb(0,0,1,0.5))
points(df[ref[ref[,i] == 2,]$FID,paste0(i,'.contrast')],
       df[ref[ref[,i] == 2,]$FID,paste0(i,'.size')],
       cex=0.5, pch=16,col=rgb(1,0,0,0.5))
points(df[ref[is.na(ref[,i]),]$FID,paste0(i,'.contrast')],
       df[ref[is.na(ref[,i]),]$FID,paste0(i,'.size')],
       cex=0.7, pch="x",col=rgb(0,0,0,0.5))
legend('bottomright',legend=c(paste0("AA (",table(ref[,i])[3],")"),paste0("AB (",table(ref[,i])[2],")"),paste0("BB (",table(ref[,i])[1],")"),paste0("NA (",sum(is.na(ref[,i])),")")),
       col = c("red","blue","green","grey"),
       #lty = c(16,16,16,NA),
       pch = c(16,16,16,4),cex = 1,
       title = "Allele(# of sample)"
       )
?legend()
dev.off()
ref[ref[,i] == 0,]$FID
i = 'AX.106718326.contrast'
i = str_sub(i,1,-10)
i
summary()
table(ref$AX.11373064)
table(ref$AX.12588615)[3]
sum(is.na(ref$AX.12588615))
setwd("")
pdf()
png("plot/test.pnd",height = 600,width = 600)
for (i in 2:ncol(df)) {
  
}
head(ref)[1:10]
target = colnames(ref)[i]
for (i in 7:ncol(ref)) {
#for (i in 7:10) {
  target = colnames(ref)[i]
  png(paste0("plot/",target,".png"),height = 600,width = 600)
  plot(df[,paste0(target,".contrast")],df[,paste0(target,'.size')], cex=0.5, pch=0,col=rgb(1,1,1,0.5),
       xlab = "Contrast = log2(SNP-A / SNP-B)",ylab = "Size = (log2(SNP-A) + log2(SNP-B))/2",
      main= paste0("ID : ",target))
  #dev.off()
  #plot(xlab = "Contrast = log2(SNP-A / SNP-B)",ylab = "Size = (log2(SNP-A) + log2(SNP-B))/2")
  points(df[ref[ref[,target] == 0,]$FID,paste0(target,'.contrast')],
         df[ref[ref[,target] == 0,]$FID,paste0(target,'.size')],
         cex=0.5, pch=16,col=rgb(0,1,0,0.5))
  points(df[ref[ref[,target] == 1,]$FID,paste0(target,'.contrast')],
         df[ref[ref[,target] == 1,]$FID,paste0(target,'.size')],
         cex=0.5, pch=16,col=rgb(0,0,1,0.5))
  points(df[ref[ref[,target] == 2,]$FID,paste0(target,'.contrast')],
         df[ref[ref[,target] == 2,]$FID,paste0(target,'.size')],
         cex=0.5, pch=16,col=rgb(1,0,0,0.5))
  points(df[ref[is.na(ref[,target]),]$FID,paste0(target,'.contrast')],
         df[ref[is.na(ref[,target]),]$FID,paste0(target,'.size')],
         cex=0.7, pch="x",col=rgb(0,0,0,0.5))
  legend('bottomright',legend=c(paste0("AA (",table(ref[,target])[1],")"),paste0("AB (",table(ref[,target])[2],")"),paste0("BB (",table(ref[,target])[3],")"),paste0("NA (",sum(is.na(ref[,target])),")")),
         col = c("red","blue","green","grey"),
         #lty = c(16,16,16,NA),
         pch = c(16,16,16,4),cex = 1,
         title = "Allele(# of sample)"
        )
  
  dev.off()
  
}

#### 회사별

DNAlink <- read.table("~/Desktop/KCDC/KKY/00.sampleInfo/DANlink.2020.cel.list.txt")
tera <- read.table("~/Desktop/KCDC/KKY/00.sampleInfo/2020.7th.tera.cel.list.txt")
DNAlink$company <- "DNAlink"
tera$company <- "Teragen"
rmPCA <- read.table("~/Desktop/KCDC/KKY/01.1stQC.PCA/ALL/notchr6_14/rmPCA.txt")
rmMH <- read.table("~/Desktop/KCDC/KKY/01.1stQC.PCA/ALL/notchr6_14/rmLQSamples.txt")
head(rmPCA)
head(rmMH)

#ref2 <- rbind(DNAlink,tera)
head(df[1:5])
head(ref2)
target <- 'AX.106718326'
plot(df[,'AX.106718326.contrast'],df[,'AX.106718326.size'], cex=0.5, pch=0,col=rgb(1,1,1,0.5),
     xlab = "Contrast = log2(SNP-A / SNP-B)",ylab = "Size = (log2(SNP-A) + log2(SNP-B))/2",
)
points(df[df$CEL %in% DNAlink$V1,paste0(target,'.contrast')],
       df[df$CEL %in% DNAlink$V1,paste0(target,'.size')],
       cex=0.5, pch=16,col=rgb(0,0,1,0.5))
points(df[df$CEL %in% tera$V1,paste0(target,'.contrast')],
       df[df$CEL %in% tera$V1,paste0(target,'.size')],
       cex=0.5, pch=16,col=rgb(1,0,0,0.5))
points(df[df$ID %in% rmPCA$V1,paste0(target,'.contrast')],
       df[df$ID %in% rmPCA$V1,paste0(target,'.size')],
       cex=0.5, pch=16,col=rgb(0,1,0,1))
points(df[df$ID %in% rmMH$V1,paste0(target,'.contrast')],
       df[df$ID %in% rmMH$V1,paste0(target,'.size')],
       cex=0.5, pch=16,col=rgb(0,0,0,1))


legend('bottomright',legend=c("Teragen","DNAlink","missing-het(44)","PCA(159)"),
       col = c("red","blue","black","green"),
       #lty = c(16,16,16,NA),
       pch = c(16,16,16,16),cex = 1,
       #title = "Allele(# of sample)"
)

for (i in 7:ncol(ref)) {
#for (i in 7:10) {
  target = colnames(ref)[i]
  png(paste0("plot/",target,".png"),height = 600,width = 600)
  plot(df[,paste0(target,".contrast")],df[,paste0(target,'.size')], cex=0.5, pch=0,col=rgb(1,1,1,0.5),
       xlab = "Contrast = log2(SNP-A / SNP-B)",ylab = "Size = (log2(SNP-A) + log2(SNP-B))/2",
       main= paste0("ID : ",target))
  #dev.off()
  #plot(xlab = "Contrast = log2(SNP-A / SNP-B)",ylab = "Size = (log2(SNP-A) + log2(SNP-B))/2")
  points(df[df$CEL %in% DNAlink$V1,paste0(target,'.contrast')],
         df[df$CEL %in% DNAlink$V1,paste0(target,'.size')],
         cex=1, pch=16,col=rgb(0,0,1,0.5))
  points(df[df$CEL %in% tera$V1,paste0(target,'.contrast')],
         df[df$CEL %in% tera$V1,paste0(target,'.size')],
         cex=1, pch=16,col=rgb(1,0,0,0.5))
  points(df[df$ID %in% rmPCA$V1,paste0(target,'.contrast')],
         df[df$ID %in% rmPCA$V1,paste0(target,'.size')],
         cex=0.7, pch=16,col=rgb(0,1,0,1))
  points(df[df$ID %in% rmMH$V1,paste0(target,'.contrast')],
         df[df$ID %in% rmMH$V1,paste0(target,'.size')],
         cex=0.7, pch=16,col=rgb(0,0,0,1))
  
  
  legend('bottomright',legend=c("Teragen","DNAlink","missing-het(44)","PCA(159)"),
         col = c("red","blue","black","green"),
         #lty = c(16,16,16,NA),
         pch = c(16,16,16,16),cex = 1,
         #title = "Allele(# of sample)"
  )
  dev.off()
  
}
