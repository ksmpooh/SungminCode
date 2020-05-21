getwd()
setwd("C:/Users/user/Desktop/KCDC/Gastric/")


case <- read.table("plotDATA/CASE.frq",header=T)
control <- read.table("plotDATA/CONTROL.frq",header=T)
data <- merge(control,case,by="SNP")
#pdf("plot/control&case_frequency.pdf",height = 10,width=10)
png("plot/control&case_frequency.png")


plot(data$MAF.x,data$MAF.y,xlab = "Control",ylab = "Case",main = "Control & Case Frequency",col = "blue")
abline(a = 0.05,b=1,col = "red",lty = 2)
abline(a = -0.05,b=1, col = "red",lty = 2)
#head(data)
points(data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.x,
       data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.y,
       col = "red",cex = 1,pch = 1)

dev.off()


rmlist <- data[data$MAF.x-data$MAF.y >= 0.05 | data$MAF.x-data$MAF.y<=-0.05,]
dim(rmlist)
