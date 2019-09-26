getwd()
setwd("C:/Users/user/Desktop/KCDC/Gastric/")
case <- read.table("CASE_ref.frq",header = T)
control <- read.table("CONTROL_ref.frq",header = T)

data <- merge(control,case, by = "SNP")
plot(data$MAF.x,data$MAF.y,xlab = "Control",ylab = "Case",main = "Control & Case Frequency",col = "blue")
abline(a = 0.05,b=1,col = "red",lty = 2)
abline(a = -0.05,b=1, col = "red",lty = 2)
#head(data)
points(data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.x,
         data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]$MAF.y,
         col = "red",cex = 1,pch = 1)

#points(data[data$MAF.y>=0.05,],col = "blue",cex =1.5,pch = 16)
?points
rm <- data[data$MAF.x-data$MAF.y >= 0.05|data$MAF.x-data$MAF.y<=-0.05,]


