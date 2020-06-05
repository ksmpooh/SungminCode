setwd("c:/Users/user/Desktop/KCDC/transplantation/re2ndQC/freq/")

case <- read.table("KR/CASE.frq",header=T)
control <- read.table("KR/CONTROL.frq",header=T)
data <- merge(control,case,by="SNP")
#pdf("plot/control&case_frequency.pdf",height = 10,width=10)
png("KR/KR.control&case_frequency.png")


plot(data$MAF.x,data$MAF.y,xlab = "Control",ylab = "Case",main = "KR.Control & Case Frequency",col = "blue")
abline(a = 0.05,b=1,col = "red",lty = 2)
abline(a = -0.05,b=1, col = "red",lty = 2)
#head(data)
points(data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.x,
       data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.y,
       col = "red",cex = 1,pch = 1)

dev.off()
rmlist <- data[data$MAF.x-data$MAF.y >= 0.05 | data$MAF.x-data$MAF.y<=-0.05,]
dim(rmlist)
head(rmlist)

write.table(rmlist$SNP,"KR/KR.freq.SNP.remove.list.txt",col.names = F,quote = F,row.names = F)

case <- read.table("LR/CASE.frq",header=T)
control <- read.table("LR/CONTROL.frq",header=T)
data <- merge(control,case,by="SNP")
#pdf("plot/control&case_frequency.pdf",height = 10,width=10)
png("LR/LR.control&case_frequency.png")
plot(data$MAF.x,data$MAF.y,xlab = "Control",ylab = "Case",main = "LR.Control & Case Frequency",col = "blue")
abline(a = 0.05,b=1,col = "red",lty = 2)
abline(a = -0.05,b=1, col = "red",lty = 2)
#head(data)
points(data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.x,
       data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.y,
       col = "red",cex = 1,pch = 1)

dev.off()

rmlist <- data[data$MAF.x-data$MAF.y >= 0.05 | data$MAF.x-data$MAF.y<=-0.05,]
dim(rmlist)
write.table(rmlist$SNP,"LR/LR.freq.SNP.remove.list.txt",col.names = F,quote = F,row.names = F)
