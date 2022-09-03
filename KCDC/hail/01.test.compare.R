#########last cbind
setwd("~/Desktop/KCDC/hail/")
d1 <- read.table("chr22_test/hail/gaws.output.T2D.tsv",header = T)
d2 <- read.table("chr22_test/epacts/test.chr22.merge.list.txt",header = T)
head(d1)
d1 <- d1[,c(1,3,4,5)]
colnames(d1) <- c("id","hail.BETA","hail.CHISQ","hail.PVALUE")


head(d2)
d2 <- d2[,c(10,12,9)]
colnames(d2) <- c("epacts.BETA","epacts.CHISQ","epacts.PVALUE")
d <- cbind(d1,d2)
head(d)
d$check.pvalue <- 0
d[d$epacts.PVALUE < 0.05,]$check.pvalue <- 1
#d[d$epacts.PVALUE >= 0.05,]$check.pvalue <- 0
table(d$check.pvalue)
plot(main = "Beta scatter plot",d$hail.BETA,d$epacts.BETA,xlab = "Hail",ylab = "epacts",col = "darkgray")
points(d[d$check.pvalue == 1,]$hail.BETA,
       d[d$check.pvalue == 1,]$epacts.BETA,
       col = "red")

points(d[d$epacts.PVALUE < 0.05,]$hail.BETA,
       d[d$epacts.PVALUE < 0.05,]$epacts.BETA,
       col = "red")
cor(d$hail.BETA,d$epacts.BETA)
d$abs <- abs(d$hail.BETA - d$epacts.BETA)
summary(d$abs)
head(d)
par(mfrow=c(1,2))
hist(main = "Y axis : 0~",d$abs,col = adjustcolor("Red",alpha = 0.3),breaks = 10,xlab = "abs(epacts.beta - hail.beta)")
hist(main = "Y axis : 0~100",d$abs,col = adjustcolor("Red",alpha = 0.3),breaks = 30,ylim = c(0,100),xlab = "abs(epacts.beta - hail.beta)")
dev.off()
plot(x = -log10(d$hail.PVALUE),y = -log10(d$epacts.PVALUE),col = adjustcolor("blue",alpha = 0.5),main = "-log10(p-value)",xlab= "hail -log10(p-value)",ylab = "epacts -log10(p-value)")
points(d[d$check.pvalue == 1,]$hail.PVALUE,
       d[d$check.pvalue == 1,]$epacts.PVALUE,
       col = "red")
points(d[d$check.pvalue != 1,]$hail.PVALUE,
       d[d$check.pvalue != 1,]$epacts.PVALUE,
       col = "darkgrey")
summary(d)
cor(d$hail.PVALUE,d$epacts.PVALUE)
##########
setwd("~/Desktop/KCDC/hail/")
library(stringr)
#d1 <- read.table("gaws.output.tsv",header = T)
d1 <- read.table("chr22_test/hail/gaws.output.T2D.tsv",header = T)
head(d1)

d1 <- d1[,c(1,2,3,4,5)]
head(str_split_fixed(d1$alleles,'"',4))[,2]
head(str_split_fixed(str_split_fixed(d1$alleles,"\"",4)[,4],"\"",2))[,1]
d1$id <- paste0(d1$locus,"_",str_split_fixed(d1$alleles,'"',4)[,2],"/",str_split_fixed(str_split_fixed(d1$alleles,"\"",4)[,4],"\"",2)[,1])
d1 <- d1[,c(6,3,4,5)]
head(d1)

#d2 <- read.table("JG.KR.ESRD.association.sub_Total.b.firth.chr12.132000001_133841815.epacts",header = T)
d2 <- read.table("chr22_test/epacts/test.chr22.merge.list.txt",header = T)
head(d2)
#head(str_split_fixed(d2$MARKER_ID,"_",4)[,1])
d2$id <- paste0(str_split_fixed(d2$MARKER_ID,"_",4)[,1],"_",str_split_fixed(d2$MARKER_ID,"_",4)[,2])

d2 <- d2[,c("id","BETA","CHISQ","PVALUE")]
colnames(d1) <- colnames(d2)
head(d1)
head(d2)
tail(d1)
tail(d2)
library(stringr)


colnames(d1) <- c("id","hail.BETA","hail.CHISQ","hail.PVALUE")

colnames(d2) <- c("id","epacts.BETA","epacts.CHISQ","epacts.PVALUE")

table(d1$id %in% d2$id)
a <-rbind()
head(a)
head(d)
d <- merge(d1,d2,by='id')
length(unique(d$id))
head(d)
plot(main = "Beta scatter plot",d$hail.BETA,d$epacts.BETA,xlab = "Hail",ylab = "epacts")
points(d[d$hail.PVALUE < 0.05 & d$epacts.PVALUE < 0.05,],col = "red",cex = 1, pch = 1)
cor.test(d$hail.BETA,d$epacts.BETA)

#write.table(d,"chr22_test/merge.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#head(abs(d$epacts.PVALUE - d$epacts.PVALUE))
d$abs <- abs(d$hail.BETA - d$epacts.BETA)

summary(d$abs)
head(d)
par(mfrow=c(1,2))

hist(main = "Y axis : 0~",d$abs,col = adjustcolor("Red",alpha = 0.3),breaks = 10,xlab = "abs(epacts.beta - hail.beta)")
hist(main = "Y axis : 0~100",d$abs,col = adjustcolor("Red",alpha = 0.3),breaks = 30,ylim = c(0,100),xlab = "abs(epacts.beta - hail.beta)")

#hist(main = "Beta histogram : abs(epacts-hail)",d$abs,col = adjustcolor("Red",alpha = 0.3),breaks = 10)
#hist(main = "Beta histogram : abs(epacts-hail)",d$abs,col = adjustcolor("Red",alpha = 0.3),breaks = 30,ylim = c(0,100))

plot(-log10(d$hail.PVALUE),-log10(d$epacts.PVALUE),col = adjustcolor("blue",alpha = 0.5),main = "-log10(p-value)",xlab= "hail p-value",ylab = "epacts p-value")

dev.off()


cor(d$hail.PVALUE,d$epacts.PVALUE)


#####################
setwd("~/Desktop/KCDC/hail/")

d1 <- read.table("gaws.output.tsv",header = T)
head(d1)
d1 <- d1[,c(1,3,4,5)]
summary(d1)
d2 <- read.table("JG.KR.ESRD.association.sub_Total.b.firth.chr12.132000001_133841815.epacts",header = T)
head(d2)
d2$locus <- paste0(d2$CHROM,":",d2$BEGIN)
d2 <- d2[,c("locus","BETA","CHISQ","PVALUE")]
colnames(d1) <- colnames(d2)
head(d1)
head(d2)
library(stringr)

d1$POS <- str_split_fixed(d1$locus,":",2)[,2]
head(d1)
d2$POS <- str_split_fixed(d2$locus,":",2)[,2]
head(d2)

plot(main = "Beta scatter plot",d1$BETA,d2$BETA,xlab = "Hail",ylab = "epacts")

par(mfrow=c(1,2))
hist(main = "Beta histogram - hail",d1$BETA,ylim = c(0,2500),col = adjustcolor("Red",alpha = 0.3))
hist(main = "Beta histogram - epacts",d2$BETA,ylim = c(0,2500),col = adjustcolor("Blue",alpha = 0.3))

dev.off()
summary(d1$BETA)
sd(d1$BETA)

summary(d2$BETA)
sd(d2$BETA)

plot(log10(d1$PVALUE),log10(d2$PVALUE),col = adjustcolor("blue",alpha = 0.5),main = "log10(p-value)",xlab= "hail p-value",ylab = "epacts p-value")
