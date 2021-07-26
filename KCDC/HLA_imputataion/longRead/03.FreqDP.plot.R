#### freqeuncy check

# type 1 : variant call
# type 2 : multi call

setwd("~/Desktop/KCDC/long_read/analysis/03.check/02.SNP/inter/freq/")


## AF

df1 <- read.table("type1.AF.txt")
df2 <- read.table("type2.AF.txt")
head(df1)
head(df2)
out <- merge(df1,df2,by='V1')
#out <- cbind(df1,df2)
plot(x = out$V2.x,y = out$V2.y,main="Alternative Allele frequency Compare",ylab = "Multi-sample Calling AF",xlab = "Variant Calling AF")
abline(a=0,b=1,col='red',lwd=1,lty = 3)  

head(out)
out$abs = abs(out$V2.x - out$V2.y)
hist(out$abs,main="abs(Variant.Calling.AF - multi.Calling.AF)")


## DP

setwd("~/Desktop/KCDC/long_read/analysis/03.check/02.SNP/inter/DP/")
df1 <- read.table("type1.DP.txt")
df2 <- read.table("type2.DP.txt")

head(df1)
head(df2)

df1[df1 == '.'] <- 0
df2[df2 == '.'] <- 0

for (i in 2:13) {
  df1[,i] <- as.integer(df1[,i])
  df2[,i] <- as.integer(df2[,i])
}

df1$sum <- rowSums(df1[,2:13])
df2$sum <- rowSums(df2[,2:13])

out <- merge(df1[,c("V1","sum")],df2[,c("V1","sum")],by = 'V1')
head(out)


plot(x = out$sum.x,y = out$sum.y,main="DP compare",ylab = "Multi-sample Calling sum(DP)",xlab = "Variant Calling sum(DP)")

